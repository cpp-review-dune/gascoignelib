#include  "stdmultilevelsolver.h"
#include  "stdtimesolver.h"
#include  "newton.h"
#include  "compose_name.h"
#include  "gascoignemultigridmesh.h"
#include  "cg.h"
#include  <iomanip>
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"

using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------------------*/

StdMultiLevelSolver::~StdMultiLevelSolver()
{
  cout << "\n************************************************************************\n\n";
  cout << "StdMultiLevelSolver\tTIME\n";
  cout << "  Residual\t\t" << _clock_residual.read() << endl;
  cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
  cout << "MATRIX\t\t\tTIME\n";

  double vm=0.;
  double il=0.;
  double so=0.;
  double ca=0., ci=0., cs=0.;
  double re=0;
  for(int level=0;level<_SP.size();level++)
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(level));
      assert(S);
      vm += S->clock_vmult();
      il += S->clock_ilu();
      so += S->clock_solve();
      ca += S->clock_computematrix();
      ci += S->clock_computeilu();
      cs += S->clock_computesolver();
      re += S->clock_residual();
    }

  cout << "  vmult\t\t\t" << vm << endl;
  cout << "  ilu\t\t\t" << il << endl;
  cout << "  solve\t\t\t" << so << endl;
  cout << "  compute matrix\t" << ca << endl;
  cout << "  compute ilu\t\t" << ci << endl;
  cout << "  compute solver\t" << cs << endl;
  cout << "VECTOR\t\t\tTIME\n";
  cout << "  residual\t\t" << re << endl;
  cout << "\n************************************************************************\n";

  //--------------------------------------------------//

  if(DataP) delete DataP; DataP=NULL;

  for(int i=0;i<_SP.size();i++) 
    { 
      if (_SP[i]) 
	{
	  delete _SP[i]; 
	  _SP[i]=NULL; 
	}
    }
   for(int i=0; i<_Interpolator.size(); i++)  
    {
      if (_Interpolator[i]) 
	{
	  delete _Interpolator[i]; 
	  _Interpolator[i]=NULL;
	}
    }
}

/*-------------------------------------------------------------*/

StdMultiLevelSolver::StdMultiLevelSolver() : 
  _paramfile(NULL), _MAP(NULL), DataP(NULL), MON(NULL),  oldnlevels(-1)
{
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::BasicInit(const MeshAgentInterface* MAP, const ParamFile* paramfile)
{
  _MAP = MAP;

  _paramfile = paramfile;
  assert(DataP==0);
  DataP = new MultiLevelSolverData(_paramfile);

  _cor.SetName("cor");
  _res.SetName("res");
  _mg0.SetName("mg0");
  _mg1.SetName("mg1");
  _cor.SetMultiLevelSolver(this);
  _res.SetMultiLevelSolver(this);
  _mg0.SetMultiLevelSolver(this);
  _mg1.SetMultiLevelSolver(this);
  RegisterVector(_cor);
  RegisterVector(_res);
  RegisterVector(_mg0);
  RegisterVector(_mg1);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  _PD = &PDX;
  for(int level=0; level<nlevels(); ++level)  
    {
      int solverlevel = nlevels()-1-level;
      assert(GetSolver(solverlevel)) ;
      GetSolver(solverlevel)->SetProblem(PDX);
    }
  DataP->nlinfo.control().matrixmustbebuild() = 1;
}

/*-------------------------------------------------------------*/

SolverInterface* StdMultiLevelSolver::NewSolver(int solverlevel) 
{ 
  if(DataP->solver=="instat")
    {
      return new StdTimeSolver;
    }
  else
    {
      return new StdSolver;
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::RegisterVector(MultiLevelGhostVector& g) 
{
  _MlVectors.insert(g);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::RegisterVectorAndMemory() 
{
  set<MultiLevelGhostVector>::const_iterator p;
  for(p=_MlVectors.begin(); p!=_MlVectors.end(); p++) 
    {
      RegisterVectorAndMemory(*p);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::RegisterVectorAndMemory(const MultiLevelGhostVector& g) 
{
  if(nlevels()!=_SP.size())
    {
      cerr << "StdMultiLevelSolver::RegisterVectorAndMemory(): registriere:\t"<<g<<endl;
      cerr << "Solvervector not ready\n";
      assert(0);
    }
  //     cerr << "*************registriere:\t"<<g<<endl;
  _MlVectors.insert(g);
  for(int level=0; level<nlevels(); ++level)  
    {
      set<MultiLevelGhostVector>::const_iterator p = _MlVectors.begin();
      while(p!=_MlVectors.end()) 
	{
	  GetSolver(level)->RegisterVector(*p++);
	}
    }
  ReInitVector();
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewSolvers()
{
  oldnlevels = _SP.size();

  int nl = GascoigneMath::max_int(nlevels(),oldnlevels);
  _SP.resize(nl,NULL);
  ComputeLevel = _SP.size()-1;

  for(int level=0; level<nlevels(); ++level)  
    {
      const MeshInterface* MIP = GetMeshAgent()->GetMesh(level);
      assert(MIP);

      int solverlevel = nlevels()-1-level;

      // new Solvers
      if(GetSolver(solverlevel)==NULL) 
	{
	  GetSolverPointer(solverlevel) = NewSolver(solverlevel);
 	  GetSolver(solverlevel)->BasicInit(solverlevel,_paramfile,MIP);
	}
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::RegisterMatrix()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->RegisterMatrix();
    }
}

void StdMultiLevelSolver::ReInitMatrix()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitMatrix();
    }
}

void StdMultiLevelSolver::ReInitVector()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitVector();
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolverNewMesh()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      const MeshInterface* MIP = GetMeshAgent()->GetMesh(level);
      assert(MIP);

      int solverlevel = nlevels()-1-level;
      GetSolver(solverlevel)->NewMesh(solverlevel,MIP);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewMgInterpolator()
{
  for (int i=0;i<_Interpolator.size();++i)
    {
      assert(_Interpolator[i]!=NULL);
      delete _Interpolator[i];
      _Interpolator[i]=NULL;
    }
  _Interpolator.resize(nlevels()-1,NULL);

  for(int l=0; l<nlevels()-1; ++l)  
    {
      //_Interpolator[l] = new MgInterpolatorMatrix;
      _Interpolator[l] = new MgInterpolatorNested;
    }
  //
  // Interpolator [l] :   interpoliert   (l+1)->l  (fein->grob)
  //
  for (int level=0;level<nlevels()-1;++level)
    {
      int sl = nlevels()-level-2;

      const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
      assert(MT);
      assert(_Interpolator[level]);
      GetSolver(level)->ConstructInterpolator(_Interpolator[level],MT);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::ReInit(const ProblemDescriptorInterface& PDX)
{
  DataP->countresidual = 0;
  DataP->nlinfo.control().matrixmustbebuild() = 1;

  NewSolvers();
  SolverNewMesh();
  NewMgInterpolator();
  SetProblem(PDX);
  RegisterMatrix();
  RegisterVectorAndMemory();

  ReInitMatrix();
  ReInitVector();
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::InterpolateSolution(MultiLevelGhostVector& u, const GlobalVector& uold) const
{
  if(oldnlevels<=0) return;

  GetSolver(FinestLevel())->InterpolateSolution(u(FinestLevel()),uold);

  for(int l=0; l<FinestLevel(); l++)
    {
      u.Vector(l).zero();
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::vmulteqgmres(MultiLevelGhostVector& y, const MultiLevelGhostVector& x) const
{
  GetSolver(ComputeLevel)->vmulteqgmres(y(ComputeLevel),x(ComputeLevel));
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::precondition(MultiLevelGhostVector& x, MultiLevelGhostVector& y)
{
  int clevel=GascoigneMath::max_int(DataP->coarselevel,0);
  if(DataP->coarselevel==-1) clevel = FinestLevel(); 

  DataP->precinfo.reset();
  DataP->precinfo.check(0.,0.);
//   LinearMg(FinestLevel(),clevel,x,y, DataP->precinfo);
  LinearMg(ComputeLevel,clevel,x,y, DataP->precinfo);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::LinearMg(int finelevel, int coarselevel, MultiLevelGhostVector& u, const MultiLevelGhostVector& f, CGInfo& info)
{
  int clevel=coarselevel;
  for(int level=coarselevel;level<nlevels();level++)
    {
      if(GetSolver(level)->DirectSolver()) clevel=level;
    }
  // wir haben auf einem hoeheren level einen direkten loeser...
  if(clevel>finelevel)  {clevel=finelevel;}

  assert(finelevel>=clevel);

  int nl = nlevels();
  nvector<double> res(nl,0.), rw(nl,0.);

  _mg0.Vector(finelevel).equ(1.,f.Vector(finelevel));
  
  bool reached = 0;

  for(int it=0; !reached; it++)
    {
      string p = DataP->mgtype;
      string p0 = p;
      if(p=="F") p0="W";
      mgstep(res,rw,finelevel,finelevel,clevel,p0,p,u,_mg0,_mg1);
      reached = info.check(res[finelevel],rw[finelevel]);
    }
}

/*-------------------------------------------------------------*/
 
void StdMultiLevelSolver::mgstep(vector<double>& res, vector<double>& rw, 
 int l, int finelevel, int coarselevel, string& p0, string p,
 MultiLevelGhostVector& u, MultiLevelGhostVector& b, MultiLevelGhostVector& v)
{
  if(l==coarselevel)
    {
      if(p=="F") {p0="V";}
      if(l==finelevel)
	{
	  GetSolver(l)->MatrixResidual(v, u, b);
	  res[l] = v.Vector(l).norm();
	}
      GetSolver(l)->smooth_exact(u,b,v);
    }
  else
    {
      GetSolver(l)->smooth_pre(u,b,v);
      GetSolver(l)->MatrixResidual(v,u,b);
      res[l] = v.Vector(l).norm();
      
      _Interpolator[l-1]-> restrict_zero(GetSolver(l-1)->GetGV(b),GetSolver(l)->GetGV(v));
      GetSolver(l-1)->HNDistribute(b);
      GetSolver(l-1)->SetBoundaryVectorZero(b);
      u.Vector(l-1).zero();
      
      int j = 0;
      if (p0=="V") j = 1;
      if (p0=="W") j = 2;
      if (p0=="F") j = 3;
      for (int i = 0; i<j; i++)
	{
	  mgstep(res,rw,l-1,finelevel,coarselevel,p0,p,u,b,v);
  	}
      if ((l==0)&&(p=="F")) { p0="W";}
      rw[l] = u.Vector(l-1).norm();

      v.Vector(l).zero();
      GetSolver(l-1)    -> HNAverage(u);
      _Interpolator[l-1]-> prolongate_add(GetSolver(l)->GetGV(v),GetSolver(l-1)->GetGV(u));
      GetSolver(l-1) -> HNZero(u);
      GetSolver(l)   -> HNZero(v);
	     
      GetSolver(l)   -> SetBoundaryVectorZero(v);

      u.Vector(l).add(DataP->mgomega,v.Vector(l));
    
      GetSolver(l)->smooth_post(u,b,v);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Cg(MultiLevelGhostVector& x, const MultiLevelGhostVector& f, CGInfo& info)
{
  assert(0);
//   CG<SolverInterface,StdMultiLevelSolver,GhostVector> cg(*(GetSolver(ComputeLevel)),*this);
  
//   cg.solve(x,f,info);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Gmres(MultiLevelGhostVector& x, const MultiLevelGhostVector& f, CGInfo& info)
{
  assert(0);
//   int n = DataP->gmresmemsize;

//   GMRES<SolverInterface,StdMultiLevelSolver,GhostVector> gmres(*GetSolver(ComputeLevel),*this,n);
  
//   gmres.solve(x,f,info);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonOutput(NLInfo& nlinfo) const
{
  assert(MON!=0);
  MON->nonlinear_step(nlinfo.GetLinearInfo(),nlinfo);
}

/*-------------------------------------------------------------*/
 
void StdMultiLevelSolver::NewtonVectorZero(MultiLevelGhostVector& w) const
{
  w.zero();
}

/*-------------------------------------------------------------*/
 
double StdMultiLevelSolver::NewtonResidual(MultiLevelGhostVector& y, const MultiLevelGhostVector& x,const MultiLevelGhostVector& b) const
{
  _clock_residual.start();
  DataP->countresidual++;
  y.Vector(ComputeLevel).equ(1.,b.Vector(ComputeLevel));
  GetSolver(ComputeLevel)->Residual(y(ComputeLevel),x(ComputeLevel),-1.);
  GetSolver(ComputeLevel)->SetBoundaryVectorZero(y(ComputeLevel));
  _clock_residual.stop();
  return NewtonNorm(y);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonMatrixControl(MultiLevelGhostVector& u, const NLInfo& nlinfo)
{
  MON->new_matrix() = 0;

  int nm1 = nlinfo.control().newmatrix();
  int nm2 = nlinfo.control().matrixmustbebuild();

  if (nm1+nm2==0) return;
  
  MON->new_matrix() = 1;

  AssembleMatrix(u);
  ComputeIlu(u);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonLinearSolve(MultiLevelGhostVector& x, const MultiLevelGhostVector& b, const MultiLevelGhostVector& u, CGInfo& info)
{
  assert(DataP->linearsolve=="mg");

  _clock_solve.start();
  info.reset();
  
  x.zero();

  int clevel=GascoigneMath::max_int(DataP->coarselevel,0);
  if(DataP->coarselevel==-1) clevel = FinestLevel(); 
  LinearMg(ComputeLevel,clevel,x,b, info);

  GetSolver(ComputeLevel)->PressureFilterIntegrate(x(ComputeLevel));
  _clock_solve.stop();
}

/*-------------------------------------------------------------*/

double StdMultiLevelSolver::NewtonUpdate(double& rr, MultiLevelGhostVector& x, MultiLevelGhostVector& dx, MultiLevelGhostVector& r, const MultiLevelGhostVector& f, NLInfo& nlinfo)
{
  double r0 = rr;
  const CGInfo& linfo = nlinfo.GetLinearInfo();
  bool nlok = !nlinfo.control().newmatrix();
  bool lok  = linfo.control().status()=="converged";
  bool lex  = linfo.control().status()=="exploded";

  double nn = NewtonNorm(dx);
  double nr = r.Vector(ComputeLevel).norm();

  if (nn>1.e10)  lex =1;
  if (!(nn>=0.)) lex =1;
  if (nr>1.e10)  lex =1;
  if (!(nr>=0.)) lex =1;

  if(lex)
    {
      nlinfo.control().status()="diverged";
      cerr << "linear : " << linfo.control().status() << endl;
      cerr << "nonlinear : " << nn << endl;
      return NewtonNorm(dx);
    }

  double omega = 0.7;
  double relax = 1.;
    
  string message = "";

  for(int iter=0;iter<DataP->nlinfo.user().maxrelax();iter++)
    {
      if(iter>0)
	{
	  x.Vector(ComputeLevel).add(relax*(omega-1.),dx.Vector(ComputeLevel));
	  relax *= omega;
	}
      else
	{
	  x.Vector(ComputeLevel).add(relax,dx.Vector(ComputeLevel));
	}

      NewtonResidual(r,x,f);
      rr = NewtonNorm(r);
      
      message = nlinfo.check_damping(iter,rr);
      if (message=="ok")       break;
      if (message=="continue") continue;
      if (message=="exploded") 
	{
	  x.Vector(ComputeLevel).add(-relax,dx.Vector(ComputeLevel));
	  relax = 0.;
	  cout << "Damping exploded !!!!!" << endl;
	  nlinfo.control().status() = "diverged";
	  break;
	}
    }
  return NewtonNorm(dx);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleMatrix(MultiLevelGhostVector& u)
{
  SolutionTransfer(u);
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->MatrixZero();
      GetSolver(l)->AssembleMatrix(u(l),1.);
    }
  DataP->nlinfo.control().matrixmustbebuild() = 0;
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::ComputeIlu(MultiLevelGhostVector& u)
{
  SolutionTransfer(u);
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->ComputeIlu(u(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::BoundaryInit(MultiLevelGhostVector& u) const
{
  for(int l=0;l<_SP.size();l++)
    GetSolver(l)->BoundaryInit(u(l));
}

/*-------------------------------------------------------------*/

string StdMultiLevelSolver::Solve(int level, MultiLevelGhostVector& u, const MultiLevelGhostVector& b)
{
  ComputeLevel = level;

  string status;
  if(DataP->nonlinearsolve=="newton")
    {
      GetSolver(ComputeLevel)->HNAverage(u(ComputeLevel));
      newton(*this,u,b,_res,_cor,DataP->nlinfo);
      GetSolver(ComputeLevel)->HNZero(u(ComputeLevel));
      status = DataP->nlinfo.control().status();
    }
  else
    {
      assert(0);
    }
  if (status!="converged")
    {
      DataP->nlinfo.control().matrixmustbebuild() = 1;
    }
  return status;
}

/*-------------------------------------------------------------*/

double StdMultiLevelSolver::ComputeFunctional(MultiLevelGhostVector& f, const MultiLevelGhostVector& u, const Functional* FP) const
{
  return GetSolver(ComputeLevel)->ComputeFunctional(f(ComputeLevel),u(ComputeLevel),FP);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Transfer(int high, int low, MultiLevelGhostVector& u) const
{
  for(int l=high;l>=low;l--)
    {
      Transfer(l,u.Vector(l-1),u.Vector(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(int high, int low, MultiLevelGhostVector& u) const
{
  for(int l=high;l>=low;l--)
    {
      SolutionTransfer(l,u.Vector(l-1),u.Vector(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(MultiLevelGhostVector& u) const
{
  SolutionTransfer(ComputeLevel,1,u);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const
{
  GetSolver(l)->HNAverage(uf);
  Transfer(l,ul,uf);
  GetSolver(l)->HNZero(uf);
  GetSolver(l-1)->SetBoundaryVector(ul);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const
{
  assert(_Interpolator[l-1]);
  _Interpolator[l-1]->SolutionTransfer(ul,uf);
}
