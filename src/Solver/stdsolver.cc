#include  <list>
#include  <iomanip> 

#include  "giota.h"
#include  "stdsolver.h"
#include  "timepattern.h"

#include  "pointmatrix.h"
#include  "pointilu.h"
#include  "umfilu.h"

#include  "ilupermutate.h"
#include  "cuthillmckee.h"
#include  "stopwatch.h"
#include  "gascoignevisualization.h"
#include  "backup.h"
#include  "pointfunctional.h"
#include  "visu_eps.h"

#include  "q1gls2d.h"
#include  "q1gls3d.h"

#include  "edgeinfocontainer.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

StdSolver::StdSolver() : 
  _MP(NULL), _MAP(NULL), _MIP(NULL), _PDX(NULL), _PrimalSolve(1),
  omega_domain(0.), _mylevel(-1), _directsolver(0), _ZP(NULL), _paramfile(NULL), _HM(NULL)
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

void StdSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  _PDX = &PDX;
  assert(_PDX);

  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->ncomp();

  Dat.Init(_paramfile,ncomp);

  if (_MAP==NULL)
    GetMatrixPointer() = NewMatrix(ncomp, _matrixtype);

  if (_MIP==NULL)
    GetIluPointer   () = NewIlu   (ncomp, _matrixtype);
 
  MemoryVector();
  MemoryMatrix();
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
}

/*-------------------------------------------------------*/

MeshInterpretorInterface* StdSolver::NewMeshInterpretor(int dimension, const string& discname)
{
  if(dimension==2)
    {
      if (discname=="Q1")         return new Q12d;
      else if (discname=="Q1Gls") return new Q1Gls2d;
      else    assert(0);
    }
  else if(dimension==3)
    {
      if (discname=="Q1")          return new Q13d;
      else if (discname=="Q1Gls")  return new Q1Gls3d;
      else assert(0);
    }
  else
    {
      cerr << "StdSolver::NewMeshInterpretor() dimension=\""<<dimension<<"\""<<endl;
      assert(0);
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
  if(_directsolver)             return new UmfIlu(GetMatrix());
  if(_matrixtype=="point_node") return new PointIlu(ncomp,"node");
  else                          return new PointIlu(ncomp,"component");
}

/*-------------------------------------------------------*/

void StdSolver::MemoryMatrix() 
{
  ConstructPressureFilter();
  SparseStructure SA;
  GetMeshInterpretor()->Structure(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);
};

/*-------------------------------------------------------*/

void StdSolver::MemoryVector()
{
  int ncomp = GetProblemDescriptor()->GetEquation()->ncomp();
  
  GhostVectorAgent::iterator p = _NGVA.begin();
  while(p!=_NGVA.end())
    {
      const GhostVector& gv = p->first;
//       cerr << "StdSolver::MemoryVector()\t" << gv.GetName() << endl;
      assert(gv.GetSolver()==this);
      if(p->second==NULL) 
	{
	  p->second = new CompVector<double>;
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

void StdSolver::ResizeVector(CompVector<double>* x, string type) const
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
  PressureFilterIntegrate(gu);
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
  PressureFilter(y);
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
  const Equation*        EQ = GetProblemDescriptor()->GetEquation();
  const DirichletData*   DD = GetProblemDescriptor()->GetDirichletData();
  if(DD==NULL) return;
  SetBoundaryVectorStrong(f, *BM,*EQ,*DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const Equation& EQ, const DirichletData& DD) const
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
  double omega = Dat.GetOmega();
  for(int iter=0; iter<niter; iter++)
    {
      MatrixResidual(h,x,y);
      GetIlu()->solve(h);
      x.add(omega,h);
      PressureFilterIntegrate(x);
    }
  _il.stop();
}

/*-------------------------------------------------------*/

void StdSolver::smooth_pre(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
  int niter = Dat.GetIterPre();
  smooth(niter,GetGV(x),GetGV(y),GetGV(help));
}

/*-------------------------------------------------------*/

void StdSolver::smooth_exact(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
  if(_directsolver)
    {
      _so.start();
      UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
      assert(UM);
      UM->Solve(GetGV(x),GetGV(y));
      _so.stop();
    }
  else
    {
      int niter = Dat.GetIterExact();
      smooth(niter,GetGV(x),GetGV(y),GetGV(help));
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth_post(BasicGhostVector& x, const BasicGhostVector& y, BasicGhostVector& help) const
{
  int niter = Dat.GetIterPost();
  smooth(niter,GetGV(x),GetGV(y),GetGV(help));
}

/*-------------------------------------------------------*/

void StdSolver::Residual(BasicGhostVector& gy, const BasicGhostVector& gx, double d) const
{
  Residual(GetGV(gy),GetGV(gx),d);
}

/*-------------------------------------------------------*/

void StdSolver::Residual(GlobalVector& y, const GlobalVector& x, double d) const
{
  _re.start();

  HNAverage(x);

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetMeshInterpretor()->Form(y,x,*EQ,d);

  HNZero(x);
  HNDistribute(y);
  PressureFilter(y);

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
  nvector<double> y(ncomp);
  
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
  const PointFunctional* PFP = dynamic_cast<const PointFunctional*>(FP);
  if(PFP)
    {
      return ComputePointFunctional(f,u,help,PFP);
    }
  const ResidualFunctional* RFP = dynamic_cast<const ResidualFunctional*>(FP);
  if(RFP)
    {
      return ComputeResidualFunctional(f,u,help,RFP);
    }
  assert(0);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeBoundaryFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const BoundaryFunctional* FP) const
{
  HNAverage(u);
  double J = GetMeshInterpretor()->ComputeBoundaryFunctional(u,*FP);
  HNZero(u);
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const
{
  HNAverage(u);
  double J = GetMeshInterpretor()->ComputeDomainFunctional(u,*FP);
  HNZero(u);
  //cerr << "StdSolver::ComputeDomainFunctional()\t" << J << endl;
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputePointFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const PointFunctional* FP) const
{
  f.zero();

  RhsPoint(f,FP);

  HNAverage(u);      
  double J = u*f;
  HNZero(u);
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeResidualFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const ResidualFunctional* FP) const
{
  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const DirichletData* DD = GetProblemDescriptor()->GetDirichletData();

  HNAverage(u);
  f.zero();
  Rhs(f);
  Residual(f,u,-1.);
  
  const DirichletData* ABD = FP->GetDirichletData();
  assert(ABD);

  z.zero();
  SetBoundaryVectorStrong(z,*BM,*EQ,*ABD);
  
  HNAverage(z);

  double J = z*f;
  HNZero(u);
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
  //
  // ACHTUNG !!!
  // Data muessen geaveraged (und zero) werden !!!
  //
  const RightHandSideData* RHS  = GetProblemDescriptor()->GetRightHandSideData();
  const NeumannData*       NRHS = GetProblemDescriptor()->GetNeumannData();

  if(RHS)
    {
      if(RHS->GetName()=="DiracRightHandSide")
	{
	  GetMeshInterpretor()->DiracRhs(f,*RHS,d);
	}
      else if(RHS->GetName()!="zero")
	{
	  GetMeshInterpretor()->Rhs(f,*RHS,d);
	}
    }
  if(NRHS)
    {
      if(NRHS->GetName()!="zero") 
	{
	  const Equation*         EQ   = GetProblemDescriptor()->GetEquation();
	  const BoundaryManager*  BM   = GetProblemDescriptor()->GetBoundaryManager();
	  GetMeshInterpretor()->RhsNeumann(f,*EQ,BM->GetNeumannColors(),*NRHS,d);	  
	}
    }
 HNDistribute(f);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleMatrix(BasicGhostVector& gu, double d)
{
  _ca.start();
  assert(GetMatrix());

  GlobalVector& u = GetGV(gu);
  HNAverage(gu);

  GetMeshInterpretor()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d);

  DirichletMatrix();
  HNZero(gu);

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

void StdSolver::ComputeIlu(const BasicGhostVector& gu) const
{
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
  for(int c=0;c<ncomp;c++)
    {
      double s = Dat.GetIluModify(c);
      I.modify(c,s);
    }
}

/* -------------------------------------------------------*/

void StdSolver::PermutateIlu(const CompVector<double>& u) const
{
  nvector<int> perm(GetMesh()->nnodes());
  
  iota(perm.begin(),perm.end(),0);
  
  if (Dat.GetIluSort()=="cuthillmckee")
    {
      CuthillMcKee    cmc(GetMatrix()->GetStencil());
      cmc.Permutate      (perm);
    }
  else if (Dat.GetIluSort()=="streamdirection")
    {
      StreamDirection sd (GetMesh(),GetMatrix()->GetStencil(),u);
      sd.Permutate       (perm,Dat.GetStreamDirection());
    }
  else if (Dat.GetIluSort()=="vectordirection")
    {
      VecDirection vd (GetMesh());
      vd.Permutate    (perm,Dat.GetVectorDirection());
    }
  
  if(GetIlu()) GetIlu()->ConstructStructure(perm,*GetMatrix());
}

/* -------------------------------------------------------*/

void StdSolver::Visu(const string& name, const BasicGhostVector& gu, int i) const
{
  HNAverage(gu);

  GascoigneVisualization Visu;
  Visu.SetMesh(GetMesh());  
  Visu.AddVector(&GetGV(gu));

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();

  HNZero(gu);
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

void StdSolver::ConstructPressureFilter()
{
   if(Dat.GetPfilter().size()==0) return;
   omega_domain = GetMeshInterpretor()->PressureFilter(PF);
}

/*--------------------------------------------------------*/

double StdSolver::EnergyEstimator(nvector<double>& eta, const BasicGhostVector& gu, BasicGhostVector& gf) const
{
  //GlobalVector& f = GetGV(gf);
  const GlobalVector& u = GetGV(gu);

  eta.reservesize(u.n());
  eta.zero();
  
//   const Equation*          EQ  = GetProblemDescriptor()->GetEquation();
//   const RightHandSideData* RHS = GetProblemDescriptor()->GetRightHandSideData();

  HNAverage(gu);

//   f.zero();

//   GetMeshInterpretor()->Rhs(f,*RHS,-1.);
//   GetMeshInterpretor()->Form(f,u,*EQ,1.);

  if (GetMesh()->dimension()==2)
    {
      EdgeInfoContainer<2> EIC;
      EIC.basicInit(_HM,u.ncomp());

      dynamic_cast<const Q12d*>(GetMeshInterpretor())->Jumps(EIC,u);
      dynamic_cast<const Q12d*>(GetMeshInterpretor())->JumpNorm(EIC,eta);
    }
  else if (GetMesh()->dimension()==3)
    {
      EdgeInfoContainer<3> EIC;
      EIC.basicInit(_HM,u.ncomp());

      dynamic_cast<const Q13d*>(GetMeshInterpretor())->Jumps(EIC,u);
      dynamic_cast<const Q13d*>(GetMeshInterpretor())->JumpNorm(EIC,eta);
    }
  //HNDistribute(gf);
  HNZero(gu);

  return eta.norm();
}

/*-------------------------------------------------------*/

void StdSolver::Read(BasicGhostVector& gu, const string& filename) const
{
  GlobalVector& u = GetGV(gu);
  ReadBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::Write(const BasicGhostVector& gu, const string& filename) const
{
  const GlobalVector& u = GetGV(gu);
  WriteBackUp(u,filename);
}

/*-------------------------------------------------------*/

int StdSolver::RhsPoint(GlobalVector& f, const PointFunctional* FP) const
{
  return GetMeshInterpretor()->RhsPoint(f,FP);
}

/*-------------------------------------------------------*/

void StdSolver::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  GetMeshInterpretor()->ConstructInterpolator(I,MT);
}

/* -------------------------------------------------------*/

void StdSolver::PressureFilterIntegrate(BasicGhostVector& x) const
{
  PressureFilterIntegrate(GetGV(x));
}

/* -------------------------------------------------------*/

nvector<double> StdSolver::IntegrateSolutionVector(const GlobalVector& u) const
{
  assert(PF.size());

  HNAverage(u);
  assert(GetMeshInterpretor()->HNZeroCheck(u)==0);

  nvector<double> dst(u.ncomp(),0.);
  nvector<double>::const_iterator pf = PF.begin();
  
  for (int j=0; j<u.n(); j++)
    {
      for (int c=0; c<u.ncomp(); c++)
	{
	  dst[c] += u(j,c)* *pf;
	}      
      pf++;
    }  
  HNZero(u);
  return dst;
}

/* -------------------------------------------------------*/

void StdSolver::PressureFilterIntegrate(GlobalVector& x) const
{
  // In each nonlinear step: applied to Newton correction,
  // in each smoothing step
  //
  if(Dat.GetPfilter().size()==0) return;
  HNAverage(x);

  nvector<double> mean = IntegrateSolutionVector(x);

  for (int i=0; i<Dat.GetPfilter().size(); i++)
    {
      int   comp = Dat.GetPfilter()[i];
      double sub = mean[comp]/omega_domain;
      x.CompAdd(comp,-sub);
    }  
  HNZero(x);
}

/*---------------------------------------------------*/

void StdSolver::PressureFilter(GlobalVector& x) const
{
  if (Dat.GetPfilter().size()==0) return;

  HNAverage(x);
  for (int i=0; i<Dat.GetPfilter().size(); i++)
    {
      int comp = Dat.GetPfilter()[i];
      double d = 0.;
      for (int j=0; j<x.n(); j++)
	{
	  d += x(j,comp);
	}      
      d /= GetMeshInterpretor()->n();
      
      x.CompAdd(comp,-d);
    }
  HNZero(x);
}

/*---------------------------------------------------*/

void StdSolver::setHierarchicalMeshPointer(const HierarchicalMesh* HM)
{
  _HM = HM;
}
