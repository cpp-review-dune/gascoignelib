#include  <list>
#include  <iomanip> 

#include  "giota.h"
#include  "stdsolver.h"
#include  "timepattern.h"

#include  "pointmatrix.h"
#include  "pointilu.h"
/*--------------------------------*/
#ifdef __WITH_UMFPACK__
#include  "umfilu.h"
#endif
/*--------------------------------*/
#include  "ilupermutate.h"
#include  "cuthillmckee.h"
#include  "stopwatch.h"
#include  "gascoignevisualization.h"
#include  "backup.h"
#include  "visu_eps.h"

#include  "diracrighthandside.h"

#include  "q1gls2d.h"
#include  "q1gls3d.h"

#include  "edgeinfocontainer.h"

#include  "glsequation.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
StdSolver::StdSolver() : 
  _MP(NULL), _MAP(NULL), _MIP(NULL), _PDX(NULL), _PrimalSolve(1),
  //  omega_domain(0.), 
  _mylevel(-1), _directsolver(0), _ZP(NULL), _paramfile(NULL), _HM(NULL)
{
}

/*-----------------------------------------*/

StdSolver::~StdSolver()
{
  if(_MAP) delete _MAP; _MAP=NULL;
  if(_MIP) delete _MIP; _MIP=NULL;
  if(_ZP)  delete _ZP;  _ZP=NULL;
}

/*-------------------------------------------------------*/

void StdSolver::_check_consistency(const Equation* EQ,const MeshInterpretorInterface* MP) const
{
  const GlsEquation* GLS = dynamic_cast<const GlsEquation*>(EQ);
  if(GLS)
    {
      if( (MP->GetName()!="Q1Gls2d") && (MP->GetName()!="Q1Gls3d") )
	{
	  cerr << "StdSolver::_check_consistency: not a consistent MeshInterpretor\n";
	  assert(0);
	  exit(-111);
	}
    }
}

/*-------------------------------------------------------*/

void StdSolver::MatrixZero() const
{
  GetMatrix()->zero();
}

/*-------------------------------------------------------*/

void StdSolver::OutputSettings() const
{
  cout << "=====================================" << endl;
  cout << "Solver:          " << GetName() << endl;
  cout << "MeshInterpretor: " << GetMeshInterpretor()->GetName()  << endl;
  GetProblemDescriptor()->OutputSettings(cout);
  cout << "=====================================" << endl;
}

/*-------------------------------------------------------*/

void StdSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->ncomp();
  
  if (_MAP==NULL)
    GetMatrixPointer() = NewMatrix(ncomp, _matrixtype);
  
  if (_MIP==NULL)
    GetIluPointer   () = NewIlu   (ncomp, _matrixtype);
}

/*-------------------------------------------------------*/

void StdSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  _PDX = &PDX;
  assert(_PDX);
  
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  _check_consistency(EQ,GetMeshInterpretor());
//  int ncomp = EQ->ncomp();
  
//   Dat.ReInit(_paramfile,ncomp);
}

/*-------------------------------------------------------*/

void StdSolver::NewMesh(int level, const MeshInterface* mp)
{
  _MP = mp;
  assert(_MP);

  GetMeshInterpretor()->ReInit(_MP);
}

/*-------------------------------------------------------*/

void StdSolver::BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)
{
  assert(MP);
  _MP = MP;
  _paramfile = paramfile;
  _mylevel=level;

  DataFormatHandler DFH;
  DFH.insert("matrixtype" , &_matrixtype, "point_node");
  DFH.insert("ndirect"    , &_ndirect   , 100);
  DFH.insert("disc", &_discname, "Q1");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Solver");

  if(MP->nnodes()<_ndirect) 
    {
      _directsolver=1;
    }

  int dimension = MP->dimension();

  GetMeshInterpretorPointer() = NewMeshInterpretor(dimension, _discname);
  assert(_ZP);

  GetMeshInterpretor()->BasicInit(_paramfile);

  _Dat.BasicInit(_paramfile);
  _PF.SetComponents(_Dat.GetPfilter());
}

/*-------------------------------------------------------*/

MeshInterpretorInterface* StdSolver::NewMeshInterpretor(int dimension, const string& discname)
{
  if(dimension==2)
    {
      if (discname=="Q1")         return new Q12d;
      else if (discname=="Q1Gls") return new Q1Gls2d;
      else 
	{         
	  cerr << " StdSolver::NewMeshInterpretor()\tunknown discname="<<discname<<endl;
	  abort();
	}
    }
  else if(dimension==3)
    {
      if (discname=="Q1")          return new Q13d;
      else if (discname=="Q1Gls")  return new Q1Gls3d;
      else 
	{         
	  cerr << " StdSolver::NewMeshInterpretor()\tunknown discname="<<discname<<endl;
	  abort();
	}
    }
  else
    {
      cerr << "StdSolver::NewMeshInterpretor() dimension=\""<<dimension<<"\""<<endl;
      abort();
    }
}

/*-------------------------------------------------------------*/

MatrixInterface* StdSolver::NewMatrix(int ncomp, const string& matrixtype) 
{
  if(_directsolver)            return new PointMatrix(ncomp,"node");
  if(matrixtype=="point_node") return new PointMatrix(ncomp,"node");
  else                         return new PointMatrix(ncomp,"component");
}

/*-------------------------------------------------------------*/

IluInterface* StdSolver::NewIlu(int ncomp, const string& matrixtype) 
{
#ifdef __WITH_UMFPACK__
  if(_directsolver)             return new UmfIlu(GetMatrix());
#endif
  if(_matrixtype=="point_node") return new PointIlu(ncomp,"node");
  else                          return new PointIlu(ncomp,"component");
}

/*-------------------------------------------------------*/

void StdSolver::ReInitMatrix() 
{
  GetMeshInterpretor()->InitFilter(_PF);
  SparseStructure SA;
  GetMeshInterpretor()->Structure(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector()
{
  int ncomp = GetProblemDescriptor()->GetEquation()->ncomp();
  
  GhostVectorAgent::iterator p = _NGVA.begin();
  while(p!=_NGVA.end())
    {
      const GhostVector& gv = p->first;
      assert(gv.GetSolver()==this);
      if(p->second==NULL) 
	{
	  //int n = GetMeshInterpretor()->n();
	  p->second = new GlobalVector;
	  p->second->ncomp() = ncomp;
	}
      ResizeVector(p->second,gv.GetType());
      p++;
    }
}

/*-------------------------------------------------------*/

void StdSolver::Zero(BasicGhostVector& dst) const
{
  GetGV(dst).zero();
}

/*-----------------------------------------*/

double StdSolver::NewtonNorm(const BasicGhostVector& u) const
{
  return GetGV(u).norm();
}

/*-----------------------------------------*/

void StdSolver::ResizeVector(GlobalVector* x, string type) const
{
  int n = GetMeshInterpretor()->n();
//   cerr << "StdSolver::ResizeVector() n="<<n<<endl; 
  x->reservesize(n);
}

/*-----------------------------------------*/

void StdSolver::HNAverage(const GlobalVector& x) const {
  GlobalVector& v = const_cast<GlobalVector&>(x);
  GetMeshInterpretor()->HNAverage(v);
}
void StdSolver::HNZero(const GlobalVector& x) const {
  GlobalVector& v = const_cast<GlobalVector&>(x);
  GetMeshInterpretor()->HNZero(v);
}
bool StdSolver::HNZeroCheck(const GlobalVector& x) const {
  assert(GetMeshInterpretor()->n()==x.n());
  GlobalVector& v = const_cast<GlobalVector&>(x);
  return GetMeshInterpretor()->HNZeroCheck(v);
}
void StdSolver::HNDistribute(GlobalVector& x) const {
  GetMeshInterpretor()->HNDistribute(x);
}
void StdSolver::HNAverageData() const {
  GetMeshInterpretor()->HNAverageData();
}
void StdSolver::HNZeroData() const {
  GetMeshInterpretor()->HNZeroData();
}

void StdSolver::HNAverage(const BasicGhostVector& x) const {
  HNAverage(GetGV(x));
}
void StdSolver::HNZero(const BasicGhostVector& x) const {
  HNZero(GetGV(x));
}
void StdSolver::HNDistribute(BasicGhostVector& x) const {
  HNDistribute(GetGV(x));
}

/*-------------------------------------------------------*/

void StdSolver::InterpolateSolution(BasicGhostVector& gu, const GlobalVector& uold) const
{
  GlobalVector& u = GetGV(gu);

  u.zero();
  GetMeshInterpretor()->InterpolateSolution(u, uold);
  SubtractMean(gu);
}

/*-----------------------------------------*/

void StdSolver::vmult(GlobalVector& y, const GlobalVector& x, double d) const
{
  _vm.start();
  GetMatrix()->vmult(y,x,d);
  _vm.stop();
}

/*-----------------------------------------*/

void StdSolver::vmulteq(GlobalVector& y, const GlobalVector& x, double d) const
{
  _vm.start();
  y.zero();
  GetMatrix()->vmult(y,x,d);
  _vm.stop();
}

/*-----------------------------------------*/

void StdSolver::residualgmres(BasicGhostVector& gy, const BasicGhostVector& gx, const BasicGhostVector& gb) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& x = GetGV(gx);
  const GlobalVector& b = GetGV(gb);

  vmulteq(y,x,1.);
  y.sadd(-1.,1.,b);
  SetBoundaryVectorZero(gy);
}
/*-----------------------------------------*/

void StdSolver::vmulteqgmres(BasicGhostVector& gy, const BasicGhostVector& gx) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& x = GetGV(gx);
  vmulteq(y,x,1.);
}

/*-----------------------------------------*/

void StdSolver::MatrixResidual(BasicGhostVector& gy, const BasicGhostVector& gx, const BasicGhostVector& gb) const
{
  MatrixResidual(GetGV(gy),GetGV(gx),GetGV(gb));
}

/*-----------------------------------------*/

void StdSolver::MatrixResidual(GlobalVector& y, const GlobalVector& x, const GlobalVector& b) const
{
  y.equ(1.,b);
  vmult(y,x,-1.);
  SubtractMeanAlgebraic(y);
  //SubtractMean(y);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorZero(BasicGhostVector& gf) const
{
  SetBoundaryVectorZero(GetGV(gf));
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorZero(GlobalVector& f) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletComponents(col);
      GetMeshInterpretor()->StrongDirichletVectorZero(f, col, comp);
    }
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVector(BasicGhostVector& f) const
{
  SetBoundaryVector(GetGV(f));
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVector(GlobalVector& f) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const DirichletData*   DD = GetProblemDescriptor()->GetDirichletData();
  if(DD==NULL) 
  {
    if(BM->GetDirichletColors().size()!=0) 
    {
      cerr << "No DirichetData given but DirichetColors in ParamFile!" << endl;
      abort();
    }
    return;
  }
  SetBoundaryVectorStrong(f,*BM,*DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const DirichletData& DD) const
{
  IntSet PrefCol = DD.preferred_colors();
  list<int> colors(BM.GetDirichletColors().begin(), 
		   BM.GetDirichletColors().end());
  
  for(IntSet::const_iterator p=PrefCol.begin();p!=PrefCol.end();p++)
    {
      int col = *p;
      colors.remove(col);
      colors.push_back(col);
    }
  for(list<int>::const_iterator p=colors.begin();p!=colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM.GetDirichletComponents(col);
      GetMeshInterpretor()->StrongDirichletVector(f, DD, col, comp);
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth(int niter, GlobalVector& x, const GlobalVector& y, GlobalVector& h) const
{
  _il.start();
  double omega = _Dat.GetOmega();
  for(int iter=0; iter<niter; iter++)
    {
      MatrixResidual(h,x,y);
      GetIlu()->solve(h);
      x.add(omega,h);
      SubtractMean(x);
    }
  _il.stop();
}

/*-------------------------------------------------------*/

void StdSolver::smooth_pre(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
  int niter = _Dat.GetIterPre();
  smooth(niter,GetGV(x),GetGV(y),GetGV(help));
}

/*-------------------------------------------------------*/

void StdSolver::smooth_exact(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver)
    {
      _so.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
      UM->Solve(GetGV(x),GetGV(y));
      _so.stop();
    }
  else
#endif
    {
      int niter = _Dat.GetIterExact();
      smooth(niter,GetGV(x),GetGV(y),GetGV(help));
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth_post(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
  int niter = _Dat.GetIterPost();
  smooth(niter,GetGV(x),GetGV(y),GetGV(help));
}

/*-------------------------------------------------------*/

void StdSolver::Form(BasicGhostVector& gy, const BasicGhostVector& gx, double d) const
{
  Form(GetGV(gy),GetGV(gx),d);
}

/*-------------------------------------------------------*/

void StdSolver::Form(GlobalVector& y, const GlobalVector& x, double d) const
{
  _re.start();

  HNAverage(x);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetMeshInterpretor()->Form(y,x,*EQ,d);

  const RobinData* RD = GetProblemDescriptor()->GetRobinData();
  if(RD)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetMeshInterpretor()->BoundaryForm(y,x,BM->GetRobinColors(),*RD,d);
  }

  HNZero(x);
  HNZeroData();
  HNDistribute(y);
  SubtractMeanAlgebraic(y);

  _re.stop();
}

/*-------------------------------------------------------*/

void StdSolver::BoundaryInit(BasicGhostVector& Gu)  const
{
  GlobalVector& u = GetGV(Gu);
  const DirichletData* DD = GetProblemDescriptor()->GetDirichletData();
  if(DD==NULL)
    {
      u.zero();
      cerr << "StdSolver::BoundaryInit():\t No DirichletData !!\n";
      return;
    }
  
  int color = 0;
  int ncomp = GetProblemDescriptor()->GetEquation()->ncomp();
  DoubleVector y(ncomp);
  
  int dim = GetMesh()->dimension();
  for (int ind=0; ind<GetMesh()->nnodes(); ind++)
    {
      if (dim==2) (*DD)(y,GetMesh()->vertex2d(ind),color);
      else        (*DD)(y,GetMesh()->vertex3d(ind),color);
      
      for (int c=0; c<u.ncomp(); c++)
  	{
	  u(ind,c) = y[c];
	}
    }
}

/*-------------------------------------------------------*/

void StdSolver::SolutionInit(BasicGhostVector& Gu)  const
{
  GlobalVector& u = GetGV(Gu);
  const InitialCondition* u0 = GetProblemDescriptor()->GetInitialCondition();
  assert(u0);

  int ncomp = GetProblemDescriptor()->GetEquation()->ncomp();
  assert(ncomp==u.ncomp());
  
  for (int ind=0; ind<GetMesh()->nnodes(); ind++)
    {
      if (GetMesh()->dimension()==2)
	{
	  for (int c=0; c<u.ncomp(); c++)
	    {
	      u(ind,c) = (*u0)(c,GetMesh()->vertex2d(ind));
	    }
	}
      else
	{
	  for (int c=0; c<u.ncomp(); c++)
	    {
	      u(ind,c) = (*u0)(c,GetMesh()->vertex3d(ind));
	    }
	}
    }
}

/*-------------------------------------------------------*/

void StdSolver::ComputeError(const BasicGhostVector& u, GlobalVector& err) const
{
  if(GetProblemDescriptor()->GetExactSolution()==NULL) return;
  HNAverage(u);
  GetMeshInterpretor()->ComputeError(GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeFunctional(BasicGhostVector& gf, const BasicGhostVector& gu, const Functional* FP) const
{
  return StdSolver::ComputeFunctional(GetGV(gf),GetGV(gu),FP);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeFunctional(GlobalVector& f, const GlobalVector& u, const Functional* FP) const
{
  BasicGhostVector gh("mg0");
  GlobalVector& help = GetGV(gh);

  const DomainFunctional* DFP = dynamic_cast<const DomainFunctional*>(FP);
  if(DFP)
    {
      return ComputeDomainFunctional(f,u,help,DFP);
    }
  const BoundaryFunctional* BFP = dynamic_cast<const BoundaryFunctional*>(FP);
  if(BFP)
    {
      return ComputeBoundaryFunctional(f,u,help,BFP);
    }
  const ResidualFunctional* RFP = dynamic_cast<const ResidualFunctional*>(FP);
  if(RFP)
    {
      return ComputeResidualFunctional(f,u,help,RFP);
    }
  const PointFunctional* NPFP = dynamic_cast<const PointFunctional*>(FP);
  if(NPFP)
    {
      return ComputePointFunctional(f,u,help,NPFP);
    }
  cerr << "Functional must be either of type DomainFunctional, BoundaryFunctional or PointFunctional!!!" << endl;
  abort();
}

/*-------------------------------------------------------*/

double StdSolver::ComputeBoundaryFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const BoundaryFunctional* FP) const
{
  HNAverage(u);
  HNAverageData();
  double J = GetMeshInterpretor()->ComputeBoundaryFunctional(u,*FP);
  HNZero(u);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const
{
  HNAverage(u);
  HNAverageData();
  double J = GetMeshInterpretor()->ComputeDomainFunctional(u,*FP);
  HNZero(u);
  HNZeroData();
  //cerr << "StdSolver::ComputeDomainFunctional()\t" << J << endl;
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputePointFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const PointFunctional* FP) const
{
  HNAverage(u);
  HNAverageData();
  double J = GetMeshInterpretor()->ComputePointFunctional(u,*FP);
  HNZero(u);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeResidualFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const ResidualFunctional* FP) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();

  HNAverage(u);
  HNAverageData();
  f.zero();
  Rhs(f);
  Form(f,u,-1.);
  
  const DirichletData* ABD = FP->GetDirichletData();
  assert(ABD);

  z.zero();
  SetBoundaryVectorStrong(z,*BM,*ABD);
  
  HNAverage(z);

  double J = z*f;
  HNZero(u);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

void StdSolver::Rhs(BasicGhostVector& f, double d) const
{
  StdSolver::Rhs(GetGV(f),d);
}

/*-------------------------------------------------------*/

void StdSolver::Rhs(GlobalVector& f, double d) const
{
  HNAverageData();

  const RightHandSideData* RHS  = GetProblemDescriptor()->GetRightHandSideData();
  const NeumannData*       NRHS = GetProblemDescriptor()->GetNeumannData();

  if(NRHS)
    {
       const BoundaryManager*  BM   = GetProblemDescriptor()->GetBoundaryManager();
       GetMeshInterpretor()->RhsNeumann(f,BM->GetNeumannColors(),*NRHS,d);	  
    }
  
  if(RHS)
    {
       bool done=false;
       const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
       if(DRHS)
       {
         GetMeshInterpretor()->Rhs(f,*DRHS,d);
         done = true;
       }
       const DiracRightHandSide *NDRHS = dynamic_cast<const DiracRightHandSide *>(RHS);
       if(NDRHS)
       {
         GetMeshInterpretor()->DiracRhs(f,*NDRHS,d);
         done =true;
       }
       if(!done)
       {
         cerr << "RightHandSide should be either of type DomainRightHandSide or DiracRightHandSide!!!" << endl;
         abort();
       }
    }

  HNZeroData();
  HNDistribute(f);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleMatrix(const BasicGhostVector& gu, double d)
{
  _ca.start();
  assert(GetMatrix());

  const GlobalVector& u = GetGV(gu);
  HNAverage(gu);
  HNAverageData();

  GetMeshInterpretor()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d);
  
  const RobinData* RD = GetProblemDescriptor()->GetRobinData();
  if(RD)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetMeshInterpretor()->BoundaryMatrix(*GetMatrix(),u,BM->GetRobinColors(),*RD,d);
  }

  DirichletMatrix();
  HNZero(gu);
  HNZeroData();

  _ca.stop();
}

/*-------------------------------------------------------*/

void StdSolver::DirichletMatrix() const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletComponents(col);
      GetMeshInterpretor()->StrongDirichletMatrix(*GetMatrix(), col, comp);
    }
}

/* -------------------------------------------------------*/

void StdSolver::ComputeIlu() const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver)
    {
      _cs.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
//       if(PrimalSolve==0) return;
      UM->Factorize();
      _cs.stop();
    }
  else
#endif
    {
      int ncomp = GetProblemDescriptor()->GetEquation()->ncomp();
      _ci.start();
      GetIlu()->zero();
      GetIlu()->copy_entries(GetMatrix());
      modify_ilu(*GetIlu(),ncomp);
      GetIlu()->compute_ilu();
      _ci.stop();
    }
}

/* -------------------------------------------------------*/

void StdSolver::ComputeIlu(const BasicGhostVector& gu) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver)
    {
      _cs.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
//       if(PrimalSolve==0) return;
      UM->Factorize();
      _cs.stop();
    }
  else
#endif
    {
      _ci.start();
      const GlobalVector& u = GetGV(gu);
      int ncomp = u.ncomp();
      PermutateIlu(u);
      GetIlu()->zero();
      GetIlu()->copy_entries(GetMatrix());
      modify_ilu(*GetIlu(),ncomp);
      GetIlu()->compute_ilu();
      _ci.stop();
    }
}

/*-------------------------------------------------------*/

void StdSolver::modify_ilu(IluInterface& I,int ncomp) const 
{
  if(_Dat.GetIluModify().size()==0) return;
  assert(_Dat.GetIluModify().size()==ncomp);
  for(int c=0;c<ncomp;c++)
    {
      double s = _Dat.GetIluModify(c);
      I.modify(c,s);
    }
}

/* -------------------------------------------------------*/

void StdSolver::PermutateIlu(const GlobalVector& u) const
{
  IntVector perm(GetMesh()->nnodes());
  
  iota(perm.begin(),perm.end(),0);
  
  if (_Dat.GetIluSort()=="cuthillmckee")
    {
      CuthillMcKee    cmc(GetMatrix()->GetStencil());
      cmc.Permutate      (perm);
    }
  else if (_Dat.GetIluSort()=="streamdirection")
    {
      const Equation*  EQ = GetProblemDescriptor()->GetEquation();
      assert(EQ);
      int ncomp = EQ->ncomp();
      assert(_Dat.GetStreamDirection().size()<=ncomp);
      StreamDirection sd (GetMesh(),GetMatrix()->GetStencil(),u);
      sd.Permutate       (perm,_Dat.GetStreamDirection());
    }
  else if (_Dat.GetIluSort()=="vectordirection")
    {
      VecDirection vd (GetMesh());
      vd.Permutate    (perm,_Dat.GetVectorDirection());
    }
  
  if(GetIlu()) GetIlu()->ConstructStructure(perm,*GetMatrix());
}

/* -------------------------------------------------------*/

void StdSolver::Visu(const string& name, const BasicGhostVector& gu, int i) const
{
  Visu(name,GetGV(gu),i);
}

/* -------------------------------------------------------*/

void StdSolver::Visu(const string& name, const GlobalVector& u, int i) const
{
  HNAverage(u);

  GascoigneVisualization Visu;
  Visu.SetMesh(GetMesh());  
  Visu.AddVector(&u);

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();

  HNZero(u);
}

/* -------------------------------------------------------*/

void StdSolver::VisuGrid(const string& name, int i) const
{
  assert(GetMesh());
  
  if (GetMesh()->dimension()==2)
    {
      VisuEPS eps;
      //  eps.SetOption(VisuEPS::LINEWIDTH,0.1);
      //  eps.SetOption(VisuEPS::WRITE_PATCH,1);
      eps.SetMesh(*GetMesh());
      eps.WriteGrid(name,i);
    }
}

/*-----------------------------------------*/

double StdSolver::EnergyEstimator(DoubleVector& eta, const BasicGhostVector& gu, BasicGhostVector& gf) const
{
  const GlobalVector& u = GetGV(gu);

  const Equation*          EQ  = GetProblemDescriptor()->GetEquation();
  const RightHandSideData* RHS = GetProblemDescriptor()->GetRightHandSideData();

  HNAverage(gu);
  HNAverageData();

  eta.reservesize(u.n());
  eta.zero();
  
  if (GetMesh()->dimension()==2)
    {
      EdgeInfoContainer<2> EIC;
      EIC.BasicInit(GetHierarchicalMesh(),u.ncomp());

      const Q12d* DP = dynamic_cast<const Q12d*>(GetMeshInterpretor());
      assert(DP);
      if(RHS)
      {
        const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
        assert(DRHS);
        DP->EnergyEstimator(EIC,eta,u,*EQ,*DRHS);
      }
      else
      {
        DP->EnergyEstimatorZeroRhs(EIC,eta,u,*EQ);
      }
    }
  else if (GetMesh()->dimension()==3)
    {
      EdgeInfoContainer<3> EIC;
      EIC.BasicInit(GetHierarchicalMesh(),u.ncomp());
      const Q13d* DP = dynamic_cast<const Q13d*>(GetMeshInterpretor());
      assert(DP);
      if(RHS)
      {
        const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
        assert(DRHS);
        DP->EnergyEstimator(EIC,eta,u,*EQ,*DRHS);
      }
      else
      {
        DP->EnergyEstimatorZeroRhs(EIC,eta,u,*EQ);
      }
    }
  HNZero(gu);
  HNZeroData();

  return eta.norm();
}

/*-------------------------------------------------------*/

void StdSolver::Read(BasicGhostVector& gu, const string& filename) const
{
  GlobalVector& u = GetGV(gu);
  u.zero();
  ReadBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::Write(const BasicGhostVector& gu, const string& filename) const
{
  const GlobalVector& u = GetGV(gu);
  WriteBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  GetMeshInterpretor()->ConstructInterpolator(I,MT);
}

/* -------------------------------------------------------*/

DoubleVector StdSolver::IntegrateSolutionVector(const GlobalVector& u) const
{
  HNAverage(u);
  DoubleVector dst = _PF.IntegrateVector(u);
  HNZero(u);
  return dst;
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMean(BasicGhostVector& x) const
{
  SubtractMean(GetGV(x));
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMean(GlobalVector& x) const
{
  // In each nonlinear step: applied to Newton correction,
  // in each smoothing step
  //
  if (_PF.Active())
    {
      HNZeroCheck(x);
      _PF.SubtractMean(x);
      HNZero(x);
    }
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMeanAlgebraic(BasicGhostVector& x) const
{
  SubtractMeanAlgebraic(GetGV(x));
}

/*---------------------------------------------------*/

void StdSolver::SubtractMeanAlgebraic(GlobalVector& x) const
{
  // applies to residuals
  if (_PF.Active())
    {
      HNZeroCheck(x);
      _PF.SubtractMeanAlgebraic(x);
      HNZero(x);
    }
}

/*---------------------------------------------------*/

}
