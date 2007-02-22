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
  _distribute(true), _mylevel(-1), _ndirect(1000), _directsolver(0), _discname("Q1"),
  _matrixtype("point_node"), _PrimalSolve(1), _paramfile(NULL), _useUMFPACK(true)
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

void StdSolver::_check_consistency(const Equation* EQ,const DiscretizationInterface* DI) const
{
  string eq = DI->GetName();

  bool glseq = false, glsdi = false;

  if (dynamic_cast<const GlsEquation*>(EQ))
  {
    glseq = true;
  }
  if (eq=="Q1Gls2d" || eq=="Q2Gls2d" || eq=="Q1Gls3d" || eq=="Q2Gls3d")
  {
    glsdi = true;
  }

  if(glseq && !glsdi)
  {
    cerr << "Warning: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
  }
  else if(!glseq && glsdi)
  {
    cerr << "Error: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }


  bool lpseq = false, lpsdi = false;

  if(dynamic_cast<const LpsEquation*>(EQ))
  {
    lpseq = true;
  }
  if(eq=="Q1Lps2d" || eq=="Q2Lps2d" || eq=="Q1Lps3d" || eq=="Q2Lps3d")
  {
    lpsdi = true;
  }

  if(lpseq && !lpsdi)
  {
    cerr << "Warning: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
  }
  else if(!lpseq && lpsdi)
  {
    cerr << "Error: Discretization \"" << eq << "\" doesn't go with type of given Equation!" << endl;
    abort();
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
  cout << "==================================================" << endl;
  cout << "Solver:                   " << GetName() << endl;
  cout << "Discretization:           " << GetDiscretization()->GetName()  << endl;
  GetProblemDescriptor()->OutputSettings(cout);
  cout << "==================================================" << endl;
}

/*-------------------------------------------------------*/

void StdSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();

#ifdef __WITH_UMFPACK__
  if (_useUMFPACK && _MAP!=NULL)
  {
    SimpleMatrix* SM = dynamic_cast<SimpleMatrix*>(GetMatrix());
    if ((SM && !_directsolver && _matrixtype!="point_node") || (!SM && _directsolver))
    {
      delete _MAP;
      _MAP = NULL;
    }
  }

  if (_useUMFPACK && _MIP!=NULL)
  {
    UmfIlu* UM = dynamic_cast<UmfIlu*>(GetIlu());
    if ((UM && !_directsolver) || (!UM && _directsolver))
    {
      delete _MIP;
      _MIP = NULL;
    }
  }
#endif

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
    DI.ReInit(GetMesh());
    DI.SetDataContainer(GetDiscretization()->GetDataContainer());
  }
  
  GetDiscretizationPointer() = &DI;
}

/*-------------------------------------------------------*/

void StdSolver::NewMesh(int level, const MeshInterface* mp)
{
  _MP = mp;
  assert(_MP);

  if(_MP->nnodes()<_ndirect) 
    {
      _directsolver=1;
    }
  else
    {
      _directsolver=0;
    }

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

void StdSolver::BasicInit(int level, const ParamFile* paramfile, const int dimension)
{
  _paramfile = paramfile;
  _mylevel=level;

  string xxx;

  DataFormatHandler DFH;
  DFH.insert("matrixtype" , &_matrixtype);
  DFH.insert("ndirect"    , &_ndirect);
  DFH.insert("useUMFPACK", &_useUMFPACK);
  DFH.insert("discname", &_discname);
  DFH.insert("disc", &xxx, "void");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Solver");

  if(xxx!="void")
    {
      cout << "Expression 'disc' in ParamFile not longer valid !" << endl;
      abort();
    }

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
  if(_directsolver && _useUMFPACK)             return new UmfIlu(GetMatrix());
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

void StdSolver::RegisterVector(const VectorInterface& g) 
{
  _NGVA.Register(g);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector(VectorInterface& dst)
{
  int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
  ReInitVector(dst,ncomp);
}

/*-------------------------------------------------------*/

void StdSolver::ReInitVector(VectorInterface& dst, int comp)
{
  int n = GetDiscretization()->n();
  int nc = GetDiscretization()->nc();

  // VectorInterface already registered ?
  //
  GhostVectorAgent::iterator p = _NGVA.find(dst);
  if (p==_NGVA.end()) 
    {
      _NGVA.Register(dst);
      p = _NGVA.find(dst);
    }
  assert(p!=_NGVA.end());
  
  // GlobalVector already registered ?
  //
  if (p->second==NULL) 
    {
      p->second = new GlobalVector;
    }

  // resize GlobalVector
  //
  p->second->ncomp()=comp;

  if(p->first.GetType()=="node")
    {
      p->second->reservesize(n);
    }
  else if(p->first.GetType()=="cell")
    {
      p->second->reservesize(nc);
    }
  else
  {
    cerr << "No such vector type: " << p->first.GetType() << endl;
    abort();
  }
}

/*-------------------------------------------------------*/

void StdSolver::Zero(VectorInterface& dst) const
{
  GetGV(dst).zero();
}

/*-----------------------------------------*/

double StdSolver::NewtonNorm(const VectorInterface& u) const
{
  return GetGV(u).norm_l8();
}

/*-----------------------------------------*/

void StdSolver::HNAverageData() const {
  GetDiscretization()->HNAverageData();
}
void StdSolver::HNZeroData() const {
  GetDiscretization()->HNZeroData();
}
void StdSolver::HNAverage(const VectorInterface& x) const {
  GetDiscretization()->HNAverage(const_cast<GlobalVector&>(GetGV(x)));
}
void StdSolver::HNZero(const VectorInterface& x) const {
  GetDiscretization()->HNZero(const_cast<GlobalVector&>(GetGV(x)));;
}
void StdSolver::HNDistribute(VectorInterface& x) const {
  if(GetDistribute())
  {
    GetDiscretization()->HNDistribute(GetGV(x));
  }
}

/*-------------------------------------------------------*/

void StdSolver::InterpolateSolution(VectorInterface& gu, const GlobalVector& uold) const
{
  GlobalVector& u = GetGV(gu);

  u.zero();
  GetDiscretization()->InterpolateSolution(u, uold);
  SubtractMean(gu);
}

/*-----------------------------------------*/

void StdSolver::residualgmres(VectorInterface& gy, const VectorInterface& gx, const VectorInterface& gb) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& b = GetGV(gb);

  vmulteq(gy,gx,1.);
  y.sadd(-1.,1.,b);
  SetBoundaryVectorZero(gy);
}

/*-----------------------------------------*/

void StdSolver::vmult(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _vm.start();
  GetMatrix()->vmult(GetGV(gy),GetGV(gx),d);
  _vm.stop();
}

/*-----------------------------------------*/

void StdSolver::vmulteq(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _vm.start();
  Zero(gy);
  _vm.stop();
  vmult(gy,gx,d);
}

/*-----------------------------------------*/

void StdSolver::MatrixResidual(VectorInterface& gy, const VectorInterface& gx, const VectorInterface& gb) const
{
  Equ(gy,1.,gb);
  vmult(gy,gx,-1.);
  SubtractMeanAlgebraic(gy);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorZero(VectorInterface& gf) const
{
  GlobalVector& f = GetGV(gf);

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

void StdSolver::SetBoundaryVector(VectorInterface& gf) const
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
  SetBoundaryVectorStrong(gf,*BM,*DD);
}

/*-------------------------------------------------------*/

void StdSolver::SetBoundaryVectorStrong(VectorInterface& gf, const BoundaryManager& BM, const DirichletData& DD, double d) const
{
  GlobalVector& f = GetGV(gf);

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
      GetDiscretization()->StrongDirichletVector(f, DD, col, comp, d);
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const
{
  _il.start();
  double omega = _Dat.GetOmega();
  for(int iter=0; iter<niter; iter++)
    {
      if (_Dat.GetLinearSmooth()=="ilu")
	{
	  MatrixResidual(h,x,y);
	  GetIlu()->solve(GetGV(h));
	  Add(x,omega,h);
	}
      else if (_Dat.GetLinearSmooth()=="jacobi")
	{
	  MatrixResidual(h,x,y);
	  GetMatrix()->Jacobi(GetGV(h));
	  Add(x,omega,h);
	}
      else if (_Dat.GetLinearSmooth()=="richardson")
	{
	  MatrixResidual(h,x,y);
	  Add(x,omega,h);
	}
      else if (_Dat.GetLinearSmooth()=="none")
	{}
      else
	{
	  cerr << "Smoother: " << _Dat.GetLinearSmooth() << " not valid!\n";
	  abort();
	}
      SubtractMean(x);
    }
  _il.stop();
}

/*-------------------------------------------------------*/

void StdSolver::smooth_pre(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
  int niter = _Dat.GetIterPre();
  smooth(niter,x,y,help);
}

/*-------------------------------------------------------*/

void StdSolver::smooth_exact(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver&&_useUMFPACK)
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
      smooth(niter,x,y,help);
    }
}

/*-------------------------------------------------------*/

void StdSolver::smooth_post(VectorInterface& x, const VectorInterface& y, VectorInterface& help) const
{
  int niter = _Dat.GetIterPost();
  smooth(niter,x,y,help);
}

/*-------------------------------------------------------*/

void StdSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  _re.start();

  HNAverage(gx);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetDiscretization()->Form(GetGV(gy),GetGV(gx),*EQ,d);

  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if(BE)
  {
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryForm(GetGV(gy),GetGV(gx),BM->GetBoundaryEquationColors(),*BE,d);
  }

  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);


  _re.stop();
}

/*-------------------------------------------------------*/

void StdSolver::AdjointForm(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  GlobalVector&       y = GetGV(gy);
  const GlobalVector& x = GetGV(gx);
  
  HNAverage(gx);
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

  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);
}

/*-------------------------------------------------------*/

void StdSolver::BoundaryInit(VectorInterface& Gu)  const
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

void StdSolver::SolutionInit(VectorInterface& Gu)  const
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

void StdSolver::ComputeError(const VectorInterface& u, GlobalVector& err) const
{
  if(GetProblemDescriptor()->GetExactSolution()==NULL) return;
  HNAverage(u);
  GetDiscretization()->ComputeError(GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleError(GlobalVector& eta, const VectorInterface& u, GlobalVector& err) const
{
  if(GetProblemDescriptor()->GetExactSolution()==NULL) return;
  HNAverage(u);
  GetDiscretization()->AssembleError(eta,GetGV(u),err,GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

double StdSolver::ComputeFunctional(VectorInterface& gf, const VectorInterface& gu, const Functional* FP) const
{
  VectorInterface gh("mg0");

  const DomainFunctional* DFP = dynamic_cast<const DomainFunctional*>(FP);
  if(DFP)
    {
      return ComputeDomainFunctional(gf,gu,gh,DFP);
    }
  const BoundaryFunctional* BFP = dynamic_cast<const BoundaryFunctional*>(FP);
  if(BFP)
    {
      return ComputeBoundaryFunctional(gf,gu,gh,BFP);
    }
  const ResidualFunctional* RFP = dynamic_cast<const ResidualFunctional*>(FP);
  if(RFP)
    {
      return ComputeResidualFunctional(gf,gu,gh,RFP);
    }
  const PointFunctional* NPFP = dynamic_cast<const PointFunctional*>(FP);
  if(NPFP)
    {
      return ComputePointFunctional(gf,gu,gh,NPFP);
    }
  cerr << "Functional must be either of type DomainFunctional, BoundaryFunctional or PointFunctional!!!" << endl;
  abort();
}

/*-------------------------------------------------------*/

double StdSolver::ComputeBoundaryFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const BoundaryFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  double J = GetDiscretization()->ComputeBoundaryFunctional(GetGV(gu),BM->GetBoundaryFunctionalColors(),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeDomainFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const DomainFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputeDomainFunctional(GetGV(gu),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputePointFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const PointFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputePointFunctional(GetGV(gu),*FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double StdSolver::ComputeResidualFunctional(VectorInterface& gf, const VectorInterface& gu, VectorInterface& gz, const ResidualFunctional* FP) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();

  HNAverage(gu);
  HNAverageData();
  Zero(gf);
  Rhs(gf);
  Form(gf,gu,-1.);
  
  const DirichletData* ABD = FP->GetDirichletData();
  assert(ABD);

  Zero(gz);
  SetBoundaryVectorStrong(gz,*BM,*ABD);
  
  HNAverage(gz);

  double J = ScalarProduct(gz,gf);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateCellRightHandSide(VectorInterface& f, const DomainRightHandSide& CF, double d) const
{
  assert(f.GetType()=="cell");
  HNAverageData();
  
  GetDiscretization()->EvaluateCellRightHandSide(GetGV(f),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::EvaluateBoundaryCellRightHandSide(VectorInterface& f, const BoundaryRightHandSide& CF,const BoundaryManager& BM, double d) const
{
  assert(f.GetType()=="cell");
  HNAverageData();
  
  GetDiscretization()->EvaluateBoundaryCellRightHandSide(GetGV(f),BM.GetBoundaryRightHandSideColors(),CF,d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::InterpolateDomainFunction(VectorInterface&  f, const DomainFunction& DF) const
{
  HNAverageData();
  
  if(f.GetType()=="node")
  {
    GetDiscretization()->InterpolateDomainFunction(GetGV(f),DF);
  }
  else if(f.GetType()=="cell")
  {
    GetDiscretization()->InterpolateCellDomainFunction(GetGV(f),DF);
  }
  else
  {
    cerr << "No such vector type: " << f.GetType() << endl;
    abort();
  }

  HNZeroData();
}

/*-------------------------------------------------------*/

void StdSolver::Rhs(VectorInterface& gf, double d) const
{
  GlobalVector& f = GetGV(gf);
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
  HNDistribute(gf);
}

/*-------------------------------------------------------*/

void StdSolver::AssembleMatrix(const VectorInterface& gu, double d)
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
  if(_directsolver&&_useUMFPACK)
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
      IntVector perm(GetIlu()->n());
      iota(perm.begin(),perm.end(),0);
      GetIlu()->ConstructStructure(perm,*GetMatrix());
      GetIlu()->zero();
      GetIlu()->copy_entries(GetMatrix());
      int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
      modify_ilu(*GetIlu(),ncomp);
      GetIlu()->compute_ilu();
      _ci.stop();
    }
}

/* -------------------------------------------------------*/

void StdSolver::ComputeIlu(const VectorInterface& gu) const
{
#ifdef __WITH_UMFPACK__
  if(_directsolver&&_useUMFPACK)
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
      PermutateIlu(gu);
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
  if( _Dat.GetIluModify().size()!=ncomp ) {
    cerr << "ERROR: _Dat.GetIluModify().size()="<< _Dat.GetIluModify().size() << " and ";
    cerr << "ncomp="<< ncomp << endl; 
    assert(0);
    // assert(_Dat.GetIluModify().size()==ncomp);
  }

  for(int c=0;c<ncomp;c++)
    {
      double s = _Dat.GetIluModify(c);
      I.modify(c,s);
    }
}

/* -------------------------------------------------------*/

void StdSolver::PermutateIlu(const VectorInterface& gu) const
{
  const GlobalVector& u = GetGV(gu);
  int n = GetIlu()->n();
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
  GetIlu()->ConstructStructure(perm,*GetMatrix());
}

/* -------------------------------------------------------*/

void StdSolver::Visu(const string& name, const VectorInterface& gu, int i) const
{
  if(gu.GetType()=="node")
  {
    PointVisu(name,GetGV(gu),i);
  }
  else if(gu.GetType()=="cell")
  {
    CellVisu(name,GetGV(gu),i);
  }
  else
  {
    cerr << "No such vector type: " << gu.GetType() << endl;
    abort();
  }
}

/* -------------------------------------------------------*/

void StdSolver::PointVisu(const string& name, const GlobalVector& u, int i) const
{
  GetDiscretization()->HNAverage(const_cast<GlobalVector&>(u)); 
  
  GascoigneVisualization Visu;
  Visu.SetMesh(GetMesh());  

  const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
  if(CI)
  {
    Visu.AddPointVector(CI,&u);
  }
  else
  {
    Visu.AddPointVector(&u);
  }

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();

  GetDiscretization()->HNZero(const_cast<GlobalVector&>(u)); 
}

/* -------------------------------------------------------*/

void StdSolver::CellVisu(const string& name, const GlobalVector& u, int i) const
{
  GascoigneVisualization Visu;
  Visu.SetMesh(GetMesh());  

  const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
  if(CI)
  {
    Visu.AddCellVector(CI,&u);
  }
  else
  {
    Visu.AddCellVector(&u);
  }

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();
}

/* -------------------------------------------------------*/

void StdSolver::VisuGrid(const string& name, int i) const
{
  assert(GetMesh());
  
  if (GetMesh()->dimension()==2)
    {
      VisuEPS eps(_paramfile);
      //  eps.SetOption(VisuEPS::LINEWIDTH,0.1);
      if(_discname[1]=='2')
      {
        eps.SetOption(VisuEPS::WRITE_PATCH,1);
      }
      eps.SetMesh(*GetMesh());
      eps.WriteGrid(name,i);
    }
}

/*-------------------------------------------------------*/

void StdSolver::Read(VectorInterface& gu, const string& filename) const
{
  GlobalVector& u = GetGV(gu);
  u.zero();
  ReadBackUp(u,filename);
}

/*-------------------------------------------------------*/

void StdSolver::Write(const VectorInterface& gu, const string& filename) const
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

DoubleVector StdSolver::IntegrateSolutionVector(const VectorInterface& gu) const
{
  HNAverage(gu);
  DoubleVector dst = _PF.IntegrateVector(GetGV(gu));
  HNZero(gu);
  return dst;
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMean(VectorInterface& gx) const
{
  GlobalVector& x = GetGV(gx);
  // In each nonlinear step: applied to Newton correction,
  // in each smoothing step
  //
  if (_PF.Active())
    {
      GetDiscretization()->HNZeroCheck(x);
      _PF.SubtractMean(x);
      HNZero(gx);
    }
}

/* -------------------------------------------------------*/

void StdSolver::SubtractMeanAlgebraic(VectorInterface& gx) const
{
  GlobalVector& x = GetGV(gx);
  
  // applies to residuals
  if (_PF.Active())
    {
      GetDiscretization()->HNZeroCheck(x);
      _PF.SubtractMeanAlgebraic(x);
      HNZero(gx);
    }
}

/*---------------------------------------------------*/

void StdSolver::DeleteVector(VectorInterface& p) const
{
  _NGVA.Delete(p);
}

/*-----------------------------------------*/

void StdSolver::Equ(VectorInterface& dst, double s, const VectorInterface& src) const
{
  GetGV(dst).equ(s,GetGV(src));
}

/*-----------------------------------------*/

void StdSolver::Add(VectorInterface& dst, double s, const VectorInterface& src) const
{
  GetGV(dst).add(s,GetGV(src));
}

/*-----------------------------------------*/

void StdSolver::SAdd(double s1,VectorInterface& dst, double s2, const VectorInterface& src) const
{
  GetGV(dst).sadd(s1,s2,GetGV(src));
}

/*-----------------------------------------*/

double StdSolver::Norm(const VectorInterface& dst) const
{
  return GetGV(dst).norm();
}

/*-----------------------------------------*/

double StdSolver::ScalarProduct(const VectorInterface& y, const VectorInterface& x) const
{
  return GetGV(y)*GetGV(x);
}

/*---------------------------------------------------*/

void StdSolver::AssembleDualMatrix(const VectorInterface& gu, double d)
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
