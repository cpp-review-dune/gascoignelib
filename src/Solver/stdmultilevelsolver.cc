#include  "stdmultilevelsolver.h"
#include  "stdtimesolver.h"
#include  "newton.h"
#include  "compose_name.h"
#include  "gascoignemultigridmesh.h"
#include  "cg.h"
#include  <iomanip>
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gmres.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{
StdMultiLevelSolver::~StdMultiLevelSolver()
{
  //ViewProtocoll();

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

void StdMultiLevelSolver::ViewProtocoll() const
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
  cout << "  smooth\t\t\t" << il << endl;
  cout << "  solve\t\t\t" << so << endl;
  cout << "  compute matrix\t" << ca << endl;
  cout << "  compute ilu\t\t" << ci << endl;
  cout << "  compute solver\t" << cs << endl;
  cout << "VECTOR\t\t\tTIME\n";
  cout << "  residual\t\t" << re << endl;
  cout << "\n************************************************************************\n";
}

/*-------------------------------------------------------------*/

StdMultiLevelSolver::StdMultiLevelSolver() : 
  _MAP(NULL), oldnlevels(-1), _paramfile(NULL), MON(NULL), DataP(NULL)
{
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::BasicInit(const MeshAgentInterface* MAP, const ParamFile* paramfile)
{
  _MAP = MAP;

  _paramfile = paramfile;
  
  if(!DataP)
  {
    DataP = new StdMultiLevelSolverData;
  }
  DataP->BasicInit(_paramfile);

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
}

/*-------------------------------------------------------------*/

SolverInterface* StdMultiLevelSolver::NewSolver(int solverlevel) 
{ 
  if(DataP->Solver()=="instat")
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
  g.SetMultiLevelSolver(this);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::RegisterVectorsOnSolvers() 
{
  assert(nlevels()==_SP.size());
  set<MultiLevelGhostVector>::const_iterator p;
  for(p=_MlVectors.begin(); p!=_MlVectors.end(); p++) 
    {
      for (int level=0; level<nlevels(); ++level)  
	{
	  GetSolver(level)->RegisterVector(*p);
	}
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewSolvers()
{
  oldnlevels = _SP.size();

  if (oldnlevels>nlevels())
    {
      for (int l=oldnlevels-1; l>=nlevels(); l--)
        {
          delete _SP[l];
          _SP[l] = NULL;
        }
    }
  _SP.resize(nlevels(),NULL);
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
  DataP->CountResidual() = 0;
  //  DataP->GetNLInfo().control().matrixmustbebuild() = 1;

  NewSolvers();
  SolverNewMesh();
  NewMgInterpolator();
  SetProblem(PDX);
  RegisterMatrix();
  RegisterVectorsOnSolvers();
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

void StdMultiLevelSolver::vmulteq(MultiLevelGhostVector& y, const MultiLevelGhostVector& x) const
{
  GetSolver(ComputeLevel)->vmulteq(y(ComputeLevel),x(ComputeLevel));
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
  DoubleVector res(nl,0.), rw(nl,0.);

  _mg0.Vector(finelevel).equ(1.,f.Vector(finelevel));
  
  bool reached = 0;

  for(int it=0; !reached; it++)
    {
      string p = DataP->MgType();
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

      u.Vector(l).add(DataP->MgOmega(),v.Vector(l));
    
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
  DataP->CountResidual()++;
  y.Vector(ComputeLevel).equ(1.,b.Vector(ComputeLevel));
  GetSolver(ComputeLevel)->Form(y(ComputeLevel),x(ComputeLevel),-1.);
  GetSolver(ComputeLevel)->SetBoundaryVectorZero(y(ComputeLevel));
  _clock_residual.stop();
  return NewtonNorm(y);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonMatrixControl(MultiLevelGhostVector& u, NLInfo& nlinfo)
{
  MON->new_matrix() = 0;

  int nm1 = nlinfo.control().newmatrix();
  int nm2 = nlinfo.control().matrixmustbebuild();

  if (nm1+nm2==0) return;
  
  MON->new_matrix() = 1;

  AssembleMatrix(u,nlinfo);
  ComputeIlu(u);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonLinearSolve(MultiLevelGhostVector& x, const MultiLevelGhostVector& b, CGInfo& info)
{
  _clock_solve.start();
  info.reset();
  x.zero();

  if (DataP->LinearSolve()=="mg")
    {
      int clevel=Gascoigne::max_int(DataP->CoarseLevel() ,0);
      if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 
      LinearMg(ComputeLevel,clevel,x,b, info);
    }
  else if (DataP->LinearSolve()=="gmres")
    {
      Gmres(x,b,info);
    }
  else
    {
      assert(0);
    }
  GetSolver(ComputeLevel)->SubtractMean(x(ComputeLevel));
  _clock_solve.stop();
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Gmres(MultiLevelGhostVector& x, const MultiLevelGhostVector& f, CGInfo& info)
{
  int n = DataP->GmresMemSize();

  StdSolver* S = dynamic_cast<StdSolver*>(GetSolver(ComputeLevel));
  GMRES<StdSolver,StdMultiLevelSolver,MultiLevelGhostVector> gmres(*S,*this,n);

  gmres.solve(x,f,info);
}

/*-------------------------------------------------------------*/

double StdMultiLevelSolver::NewtonUpdate(double& rr, MultiLevelGhostVector& x, MultiLevelGhostVector& dx, MultiLevelGhostVector& r, const MultiLevelGhostVector& f, NLInfo& nlinfo)
{
  const CGInfo& linfo = nlinfo.GetLinearInfo();
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

  for(int iter=0;iter<nlinfo.user().maxrelax();iter++)
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

      if(message=="ok" && DataP->ShowCompResiduals() ){
        int ncomps = r.Vector(ComputeLevel).ncomp();
        cout << "      residuals in components L8: [";
        for(int i=0;i<ncomps;i++){
          double res_comp = r.Vector(ComputeLevel).CompNormL8(i);
          cout << " " << res_comp; 
          if(i<ncomps-1) cout << ","; 
        }
        cout << " ]" << endl;
      }
      
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
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleMatrix(MultiLevelGhostVector& u, NLInfo& nlinfo)
{
  AssembleMatrix(u);
  nlinfo.control().matrixmustbebuild() = 0;
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::ComputeIlu()
{
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->ComputeIlu();
    }
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

string StdMultiLevelSolver::LinearSolve(int level, MultiLevelGhostVector& u, const MultiLevelGhostVector& b, CGInfo& info)
{
  ComputeLevel = level;

  GetSolver(ComputeLevel)->HNAverage(u(ComputeLevel));
  
  info.reset();
  
  int clevel=Gascoigne::max_int(DataP->CoarseLevel() ,0);
  if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

  LinearMg(ComputeLevel,clevel,u,b,info);

  GetSolver(ComputeLevel)->SubtractMean(u(ComputeLevel));

  GetSolver(ComputeLevel)->HNZero(u(ComputeLevel));
  string status = info.control().status();

  return status;
}

/*-------------------------------------------------------------*/

string StdMultiLevelSolver::Solve(int level, MultiLevelGhostVector& u, const MultiLevelGhostVector& b, NLInfo& nlinfo)
{
  ComputeLevel = level;

  string status;
  if(DataP->NonLinearSolve() == "newton")
    {
      GetSolver(ComputeLevel)->HNAverage(u(ComputeLevel));
      newton(*this,u,b,_res,_cor,nlinfo);
      GetSolver(ComputeLevel)->HNZero(u(ComputeLevel));
      return nlinfo.CheckMatrix();
    }
  else
    {
      assert(0);
    }
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
      GetSolver(l-1)->SetBoundaryVector(u);
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
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const
{
  assert(_Interpolator[l-1]);
  _Interpolator[l-1]->SolutionTransfer(ul,uf);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleDualMatrix(MultiLevelGhostVector& u)
{
  for(int l=0; l<nlevels(); l++)
    {
      GetSolver(l)->AssembleDualMatrix(u(l),1.);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::MemoryVector(MultiLevelGhostVector& v)
{
  RegisterVector(v);
  for(int l=0; l<nlevels(); ++l)  
    {
      StdSolver* S = dynamic_cast<StdSolver*>(GetSolver(l));
      S->MemoryVector(v(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::DeleteVector(MultiLevelGhostVector& v)
{
  _MlVectors.erase(v);
    
  for(int l=0; l<nlevels(); ++l)  
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->DeleteVector(&v(l));
    }
}


/*-------------------------------------------------------------*/

void StdMultiLevelSolver::precondition(MultiLevelGhostVector& x, MultiLevelGhostVector& y)
{
  CGInfo& precinfo = DataP->GetPrecInfo();
  precinfo.reset();
  precinfo.check(0.,0.);

  int clevel=Gascoigne::max_int(DataP->CoarseLevel(),0);
  if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

  LinearMg(ComputeLevel,clevel,x,y,precinfo);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Equ(MultiLevelGhostVector& dst, double s, const MultiLevelGhostVector& src) const
{
  for(int l=0; l<FinestLevel(); l++)
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->Equ(dst(l),s,src(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Zero(MultiLevelGhostVector& dst) const
{
  for(int l=0; l<FinestLevel(); l++)
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->Zero(dst(l));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AddNodeVector(const string& name, MultiLevelGhostVector& gq)
{
  Transfer(ComputeLevel,1,gq);
  for(int l=0; l<nlevels(); l++)
    {
      const GlobalVector& q = GetSolver(l)->GetGV(gq);
      GetSolver(l)->AddNodeVector(name,&q);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::DeleteNodeVector(const string& name)
{
  for(int l=0; l<nlevels(); l++)
    {
      GetSolver(l)->DeleteNodeVector(name);
    }
}

/*-------------------------------------------------------------*/

}
