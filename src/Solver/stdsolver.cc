#include  <list>
#include  <iomanip> 

#include  "giota.h"
#include  "stdsolver.h"

#include  "pointmatrix.h"
#include  "pointilu.h"

#include  "sparseblockilu.h"
#include  "fmatrixblock.h"
#include  "cfdblock3d.h"

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

#include  "q12d.h"
#include  "q22d.h"
#include  "q1gls2d.h"
#include  "q2gls2d.h"
#include  "q1lps2d.h"
#include  "q2lps2d.h"
#include  "q13d.h"
#include  "q23d.h"
#include  "q1gls3d.h"
#include  "q1lps3d.h"
#include  "q2lps3d.h"

#include "q2lps2dwithsecond.h"
#include "q22dwithsecond.h"
#include "q2lps3dwithsecond.h"
#include "q23dwithsecond.h"

#include  "glsequation.h"
#include  "lpsequation.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
StdSolver::StdSolver() : 
  _MP(NULL), _HM(NULL), _MAP(NULL), _MIP(NULL), _ZP(NULL), _PDX(NULL), 
  _mylevel(-1), _ndirect(1000), _directsolver(0), _discname("Q1"),
  _matrixtype("point_node"), _PrimalSolve(1), _paramfile(NULL)
// , omega_domain(0.) 
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

void StdSolver::_check_consistency(const Equation* EQ,const DiscretizationInterface* MP) const
{
  bool glseq = false, glsmi = false;
  string eq = MP->GetName();
  
  if (dynamic_cast<const GlsEquation*>(EQ))
    {
      glseq = true;
    }
  if (eq=="Q1Gls2d" || eq=="Q2Gls2d" || eq=="Q1Gls3d" || eq=="Q2Gls3d")
    {
    glsmi = true;
  }
  if(!(glseq && glsmi) && !(!glseq && !glsmi))
  {
    cerr << "Discretization \"" << MP->GetName() << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }
  
//   bool lpseq = false, lpsmi = false;
  
//   if(dynamic_cast<const LpsEquation*>(EQ))
//   {
//     lpseq = true;
//   }
//   if(MP->GetName()=="Q1Lps2d" || MP->GetName()=="Q2Lps2d" || MP->GetName()=="Q1Lps3d" || MP->GetName()=="Q2Lps3d")
//   {
//     lpsmi = true;
//   }
//   if(!(lpseq && lpsmi) && !(!lpseq && !lpsmi))
//   {
//     cerr << "Discretization \"" << MP->GetName() << "\" doesn't go with type of given Equation!" << endl;
//     abort();
//   }
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
  cout << "Solver:           " << GetName() << endl;
  cout << "Discretization:   " << GetDiscretization()->GetName()  << endl;
  GetProblemDescriptor()->OutputSettings(cout);
  cout << "=====================================" << endl;
}

/*-------------------------------------------------------*/

void StdSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();
  
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

  if (EQ) _check_consistency(EQ,GetDiscretization());
}

/*-------------------------------------------------------*/

void StdSolver::SetDiscretization(DiscretizationInterface& DI, bool init)
{
  if(init)
  {
    DI.BasicInit(GetParamfile());
    DI.ReInit(GetMesh());
    DI.SetGlobalData(GetDiscretization()->GetGlobalData());
  }
  
  GetDiscretizationPointer() = &DI;
}

/*-------------------------------------------------------*/

void StdSolver::NewMesh(int level, const MeshInterface* mp)
{
  _MP = mp;
  assert(_MP);

  GetDiscretization()->ReInit(_MP);
}

/*-----------------------------------------*/

void StdSolver::SetDefaultValues(string discname, string matrixtype, int ndirect)
{
  _discname   = discname;
  _matrixtype = matrixtype;
  _ndirect    = ndirect;
}


/*-------------------------------------------------------*/

void StdSolver::BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)
{
  assert(MP);
  _MP = MP;
  _paramfile = paramfile;
  _mylevel=level;

  DataFormatHandler DFH;
  DFH.insert("matrixtype" , &_matrixtype);
  DFH.insert("ndirect"    , &_ndirect);
  DFH.insert("discname", &_discname);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Solver");

  if(MP->nnodes()<_ndirect) 
    {
      _directsolver=1;
    }

  int dimension = MP->dimension();

  GetDiscretizationPointer() = NewDiscretization(dimension, _discname);
  assert(_ZP);

  GetDiscretization()->BasicInit(_paramfile);

  _Dat.BasicInit(_paramfile);
  _PF.SetComponents(_Dat.GetPfilter());
}

/*-------------------------------------------------------*/

DiscretizationInterface* StdSolver::NewDiscretization(int dimension, const string& discname)
{
  if (dimension==2)
    {
      if      (discname=="Q1")               return new Q12d;
      else if (discname=="Q2")               return new Q22d;
      else if (discname=="Q1Gls")            return new Q1Gls2d;
      else if (discname=="Q2Gls")            return new Q2Gls2d;
      else if (discname=="Q1Lps")            return new Q1Lps2d;
      else if (discname=="Q2Lps")            return new Q2Lps2d;
      else if (discname=="Q2WithSecond")     return new Q22dWithSecond;
      else if (discname=="Q2LpsWithSecond")  return new Q2Lps2dWithSecond;
    else 
        {         
          cerr << " Solver::NewDiscretization()\tunknown discname=" << discname << endl;
          abort();
        }
    }
  else if (dimension==3)
    {
      if      (discname=="Q1")               return new Q13d;
      else if (discname=="Q2")               return new Q23d;
      else if (discname=="Q1Gls")            return new Q1Gls3d;
      else if (discname=="Q1Lps")            return new Q1Lps3d;
      else if (discname=="Q2Lps")            return new Q2Lps3d;
      else if (discname=="Q2WithSecond")     return new Q23dWithSecond;
      else if (discname=="Q2LpsWithSecond")  return new Q2Lps3dWithSecond;
      else 
        {         
          cerr << " Solver::NewDiscretization()\tunknown discname=" << discname << endl;
          abort();
        }
    }
  else
    {
      cerr << " Solver::NewDiscretization()\tdimension must either 2 or 3" << endl;
      abort();
    }
}

/*-------------------------------------------------------------*/

MatrixInterface* StdSolver::NewMatrix(int ncomp, const string& matrixtype) 
{
  if( _directsolver || matrixtype=="point_node")
  {
    return new PointMatrix(ncomp,"node");
  }
  else if (matrixtype=="block")
  {
    if      (ncomp==1)  return new SparseBlockMatrix<FMatrixBlock<1> >;
    else if (ncomp==2)  return new SparseBlockMatrix<FMatrixBlock<2> >;
    else if (ncomp==3)  return new SparseBlockMatrix<FMatrixBlock<3> >;
    else if (ncomp==4)  return new SparseBlockMatrix<FMatrixBlock<4> >;
    else
    {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  }
  else if (matrixtype=="component")
  {
    return new PointMatrix(ncomp,"component");
  }
  else if (matrixtype=="cfd")
  {
    if (ncomp==4)
    {
      return new SparseBlockMatrix<CFDBlock3d>;
    }
    else
    {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  }
  else
  {
    cerr << "No such matrix type \"" << matrixtype<< "\"." << endl;
    abort();
  }
}

/*-------------------------------------------------------------*/

IluInterface* StdSolver::NewIlu(int ncomp, const string& matrixtype) 
{
#ifdef __WITH_UMFPACK__
  if(_directsolver)             return new UmfIlu(GetMatrix());
#endif
  // analog zu NewMatrix muss hier auch _directsolver eingehen, 
  // sonst gibts aerger nachher beim 
  // GetIlu()->copy_entries(GetMatrix());
  if(_directsolver || matrixtype=="point_node")  return new PointIlu(ncomp,"node");
  
  else if (matrixtype=="block")
  {
    if      (ncomp==1)  return new SparseBlockIlu<FMatrixBlock<1> >;
    else if (ncomp==2)  return new SparseBlockIlu<FMatrixBlock<2> >;
    else if (ncomp==3)  return new SparseBlockIlu<FMatrixBlock<3> >;
    else if (ncomp==4)  return new SparseBlockIlu<FMatrixBlock<4> >;
    else
      {
        cerr << "No SparseBlockIlu for " << ncomp << "components." << endl;
        abort();
      }
  }
  else if (matrixtype=="component") 
  {
    return new PointIlu(ncomp,"component");
  }
  else if (matrixtype=="cfd")
  {
    if (ncomp==4)  
    {
      return new SparseBlockIlu<CFDBlock3d>;
    }
  }
  cerr << "No such matrix type \"" << matrixtype << "and ncomp \"." << ncomp << endl;
  abort();
}

/*-------------------------------------------------------*/

void StdSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(_PF);
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);
}

/*-------------------------------------------------------*/

void StdSolver::RegisterVector(const BasicGhostVector& g) 
{
  _NGVA.Register(g,this);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector()
{
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
  
  for (GhostVectorAgent::iterator p=_NGVA.begin(); p!=_NGVA.end(); p++)
    {
      const GhostVector& gv = p->first;
      assert(gv.GetSolver()==this);
      if (p->second==NULL) 
        {
          p->second = new GlobalVector;
	  
	  if (gv.GetNcomp()==-1) gv.SetNcomp(ncomp);
          p->second->ncomp() = gv.GetNcomp();
        }
      ResizeVector(p->second,gv.GetType());
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
  return GetGV(u).norm_l8();
}

/*-----------------------------------------*/

void StdSolver::ResizeVector(GlobalVector* x, string type) const
{
  int n = GetDiscretization()->n();
//   cerr << "StdSolver::ResizeVector() n="<<n<<endl; 
  x->reservesize(n);
}

/*-----------------------------------------*/

void StdSolver::HNAverage(const GlobalVector& x) const {
  GlobalVector& v = const_cast<GlobalVector&>(x);
  GetDiscretization()->HNAverage(v);
}
void StdSolver::HNZero(const GlobalVector& x) const {
  GlobalVector& v = const_cast<GlobalVector&>(x);
  GetDiscretization()->HNZero(v);
}
bool StdSolver::HNZeroCheck(const GlobalVector& x) const {
  assert(GetDiscretization()->n()==x.n());
  GlobalVector& v = const_cast<GlobalVector&>(x);
  return GetDiscretization()->HNZeroCheck(v);
}
void StdSolver::HNDistribute(GlobalVector& x) const {
  GetDiscretization()->HNDistribute(x);
}
void StdSolver::HNAverageData() const {
  GetDiscretization()->HNAverageData();
}
void StdSolver::HNZeroData() const {
  GetDiscretization()->HNZeroData();
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
  GetDiscretization()->InterpolateSolution(u, uold);
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

void StdSolver::vmult(BasicGhostVector& gy, const BasicGhostVector& gx, double d) const
{
  vmult(GetGV(gy),GetGV(gx),d);
}

/*-----------------------------------------*/

void StdSolver::vmulteq(BasicGhostVector& gy, const BasicGhostVector& gx) const
{
  vmulteq(GetGV(gy),GetGV(gx),1.);
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
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletVectorZero(f, col, comp);
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
    if(BM->GetDirichletDataColors().size()!=0) 
    {
      cerr << "No DirichetData given but DirichetColors in ParamFile!" << endl;
      abort();
    }
    return;
  }
  SetBoundaryVectorStrong(f,*BM,*DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(BasicGhostVector& f, const BoundaryManager& BM, const DirichletData& DD) const
{
  SetBoundaryVectorStrong(GetGV(f),BM,DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const DirichletData& DD) const
{
  IntSet PrefCol = DD.preferred_colors();
  list<int> colors(BM.GetDirichletDataColors().begin(), 
		   BM.GetDirichletDataColors().end());
  
  for(IntSet::const_iterator p=PrefCol.begin();p!=PrefCol.end();p++)
    {
      int col = *p;
      colors.remove(col);
      colors.push_back(col);
    }
  for(list<int>::const_iterator p=colors.begin();p!=colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM.GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletVector(f, DD, col, comp);
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
  GetDiscretization()->Form(y,x,*EQ,d);

  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryForm(y,x,BM->GetBoundaryEquationColors(),*BE,d);
  }

  HNZero(x);
  HNZeroData();
  HNDistribute(y);
  SubtractMeanAlgebraic(y);

  _re.stop();
}

/*-------------------------------------------------------*/

void StdSolver::AdjointForm(GlobalVector& y, const GlobalVector& x, double d) const
{
  HNAverage(x);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetDiscretization()->AdjointForm(y,x,*EQ,d);

  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryForm(y,x,BM->GetBoundaryEquationColors(),*BE,d);
  }

  HNZero(x);
  HNZeroData();
  HNDistribute(y);
  SubtractMeanAlgebraic(y);
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
      abort();
    }
  
  int color = 0;
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
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
  const DomainInitialCondition* u0 = dynamic_cast<const DomainRightHandSide *>(GetProblemDescriptor()->GetInitialCondition());
  if(u0==NULL)
    {
      u.zero();
      return;
    }

  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
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
  GetDiscretization()->ComputeError(GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeFunctional(BasicGhostVector& gf, const BasicGhostVector& gu, const Functional* FP) const
{
  return ComputeFunctional(GetGV(gf),GetGV(gu),FP);
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
  double J = GetDiscretization()->ComputeBoundaryFunctional(u,*FP);
  HNZero(u);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const
{
  HNAverage(u);
  HNAverageData();
  double J = GetDiscretization()->ComputeDomainFunctional(u,*FP);
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
  double J = GetDiscretization()->ComputePointFunctional(u,*FP);
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
  Rhs(GetGV(f),d);
}

/*-------------------------------------------------------*/

void StdSolver::Rhs(GlobalVector& f, double d) const
{
  HNAverageData();

  const Application* RHS  = GetProblemDescriptor()->GetRightHandSide();
  const BoundaryRightHandSide* NRHS = GetProblemDescriptor()->GetBoundaryRightHandSide();

  if(RHS)
    {
       bool done=false;
       const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
       if(DRHS)
       {
         GetDiscretization()->Rhs(f,*DRHS,d);
         done = true;
       }
       const DiracRightHandSide *NDRHS = dynamic_cast<const DiracRightHandSide *>(RHS);
       if(NDRHS)
       {
         GetDiscretization()->DiracRhs(f,*NDRHS,d);
         done =true;
       }
       if(!done)
       {
         cerr << "RightHandSide should be either of type DomainRightHandSide or DiracRightHandSide!!!" << endl;
         abort();
       }
    }
  
  if(NRHS)
    {
      assert(NRHS->GetNcomp()==f.ncomp());
      const BoundaryManager*  BM   = GetProblemDescriptor()->GetBoundaryManager();
      GetDiscretization()->BoundaryRhs(f,BM->GetBoundaryRightHandSideColors(),*NRHS,d);	  
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

  GetDiscretization()->Matrix(*GetMatrix(),u,*GetProblemDescriptor()->GetEquation(),d);
  
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryMatrix(*GetMatrix(),u,BM->GetBoundaryEquationColors(),*BE,d);
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
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletMatrix(*GetMatrix(), col, comp);
    }
}

/* -------------------------------------------------------*/

void StdSolver::DirichletMatrixOnlyRow() const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntSet& Colors = BM->GetDirichletDataColors();
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& comp = BM->GetDirichletDataComponents(col);
      GetDiscretization()->StrongDirichletMatrixOnlyRow(*GetMatrix(), col, comp);
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
      int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
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
  int n = GetMatrix()->GetStencil()->n();
  IntVector perm(n);
  
  iota(perm.begin(),perm.end(),0);
  
  // if (0 && _mylevel==0 && _directsolver)
  //   {
  //     // don't do anything if we're on the lowest level and solving directly
  //   } else
  if (_Dat.GetIluSort()=="cuthillmckee")
    {
      CuthillMcKee    cmc(GetMatrix()->GetStencil());
      cmc.Permutate      (perm);
    }
  else if (_Dat.GetIluSort()=="streamdirection")
    {
      const Equation*  EQ = GetProblemDescriptor()->GetEquation();
      assert(EQ);
      int ncomp = EQ->GetNcomp();
      assert(_Dat.GetStreamDirection().size()<=ncomp);
      StreamDirection sd (GetMesh(),GetMatrix()->GetStencil(),u);
      sd.Permutate       (perm,_Dat.GetStreamDirection());
    }
  else if (_Dat.GetIluSort()=="vectordirection")
    {
      VecDirection vd (GetMesh());
      vd.Permutate    (perm,_Dat.GetVectorDirection());
    }
  if (GetIlu()) GetIlu()->ConstructStructure(perm,*GetMatrix());
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

  const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
  //Visu.AddVector(&u);
  Visu.AddVector(CI,&u);

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
  GetDiscretization()->ConstructInterpolator(I,MT);
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

void StdSolver::MemoryVector(BasicGhostVector& v)
{
  RegisterVector(v);
  GhostVectorAgent::iterator p = _NGVA.find(v);
  assert(p!=_NGVA.end());

  const GhostVector& gv = p->first;
  assert(gv.GetSolver()==this);
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
  if(p->second==NULL) 
    {
      p->second = new CompVector<double>;
      if (p->first.GetNcomp()==-1) p->first.SetNcomp(ncomp);
      p->second->ncomp() = p->first.GetNcomp();
    }
  ResizeVector(p->second,gv.GetType());
}

/*-----------------------------------------*/

void StdSolver::DeleteVector(BasicGhostVector* p) const
{
  if(p==NULL) return;

  GlobalVector* v = &GetGV(*p);
  delete v;
  _NGVA.erase(*p);
}

/*-----------------------------------------*/

void StdSolver::Equ(BasicGhostVector& dst, double s, const BasicGhostVector& src) const
{
  GetGV(dst).equ(s,GetGV(src));
}

/*-----------------------------------------*/

void StdSolver::Add(BasicGhostVector& dst, double s, const BasicGhostVector& src) const
{
  GetGV(dst).add(s,GetGV(src));
}

/*-----------------------------------------*/

double StdSolver::ScalarProduct(const BasicGhostVector& y, const BasicGhostVector& x) const
{
  return GetGV(y)*GetGV(x);
}

/*---------------------------------------------------*/

void StdSolver::AssembleDualMatrix(const BasicGhostVector& gu, double d)
{
  _ca.start();

  MatrixInterface* M = GetMatrix();

  assert(M);

  HNAverage(gu);

  const Equation& EQ = *GetProblemDescriptor()->GetEquation();
  M->zero();
  GetDiscretization()->Matrix(*M,GetGV(gu),EQ,d);
  M->transpose();

  DirichletMatrixOnlyRow();
  HNZero(gu);

  _ca.stop();
}

}
