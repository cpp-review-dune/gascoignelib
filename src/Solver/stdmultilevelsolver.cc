#include  "stdmultilevelsolver.h"
#include  "stdtimesolver.h"
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
_MAP(NULL), _cor("cor"), _res("res"), _mg0("mg0"), _mg1("mg1"),
oldnlevels(-1), _paramfile(NULL), MON(NULL), DataP(NULL)
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

void StdMultiLevelSolver::RegisterVectors() 
{
  assert(nlevels()==_SP.size());
  for (int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->RegisterVector(_cor);
      GetSolver(level)->RegisterVector(_res);
      GetSolver(level)->RegisterVector(_mg0);
      GetSolver(level)->RegisterVector(_mg1);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::ReInitVector(VectorInterface& v, int comp)
{
  for (int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitVector(v,comp);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::ReInitVector(VectorInterface& v)
{
  for (int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitVector(v);
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
          GetSolver(solverlevel)->BasicInit(solverlevel,_paramfile,GetMeshAgent()->GetDimension());
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
  RegisterVectors();
  ReInitMatrix();

  ReInitVector(_cor);
  ReInitVector(_res);
  ReInitVector(_mg0);
  ReInitVector(_mg1);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const
{
  if(oldnlevels<=0) return;

  GetSolver(FinestLevel())->InterpolateSolution(u,uold);

  for(int l=0; l<FinestLevel(); l++)
    {
      GetSolver(l)->GetGV(u).zero();
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::vmulteq(VectorInterface& y, const VectorInterface& x) const
{
  GetSolver(ComputeLevel)->vmulteq(y,x,1.);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::LinearMg(int finelevel, int coarselevel, VectorInterface& u, const VectorInterface& f, CGInfo& info)
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

  GetSolver(finelevel)->GetGV(_mg0).equ(1.,GetSolver(finelevel)->GetGV(f));
  
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
 VectorInterface& u, VectorInterface& b, VectorInterface& v)
{
  if(l==coarselevel)
    {
      if(p=="F") {p0="V";}
      if(l==finelevel)
        {
          GetSolver(l)->MatrixResidual(v, u, b);
          res[l] = GetSolver(l)->GetGV(v).norm();
        }
      GetSolver(l)->smooth_exact(u,b,v);
    }
  else
    {
      GetSolver(l)->smooth_pre(u,b,v);
      GetSolver(l)->MatrixResidual(v,u,b);
      res[l] = GetSolver(l)->GetGV(v).norm();
      
      _Interpolator[l-1]-> restrict_zero(GetSolver(l-1)->GetGV(b),GetSolver(l)->GetGV(v));
      GetSolver(l-1)->HNDistribute(b);
      GetSolver(l-1)->SetBoundaryVectorZero(b);
      GetSolver(l-1)->GetGV(u).zero();
      
      int j = 0;
      if (p0=="V") j = 1;
      if (p0=="W") j = 2;
      if (p0=="F") j = 3;
      for (int i = 0; i<j; i++)
        {
          mgstep(res,rw,l-1,finelevel,coarselevel,p0,p,u,b,v);
        }
      if ((l==0)&&(p=="F")) { p0="W";}
      rw[l] = GetSolver(l-1)->GetGV(u).norm();

      GetSolver(l)->GetGV(v).zero();
      GetSolver(l-1)    -> HNAverage(u);
      _Interpolator[l-1]-> prolongate_add(GetSolver(l)->GetGV(v),GetSolver(l-1)->GetGV(u));
      GetSolver(l-1) -> HNZero(u);
      GetSolver(l)   -> HNZero(v);
     
      GetSolver(l)   -> SetBoundaryVectorZero(v);

      GetSolver(l)->GetGV(u).add(DataP->MgOmega(),GetSolver(l)->GetGV(v));
    
      GetSolver(l)->smooth_post(u,b,v);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Cg(VectorInterface& x, const VectorInterface& f, CGInfo& info)
{
  assert(0);
//   CG<SolverInterface,StdMultiLevelSolver,GhostVector> cg(*(GetSolver(ComputeLevel)),*this);
  
//   cg.solve(x,f,info);
}


/*-------------------------------------------------------------*/

void StdMultiLevelSolver::newton(VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info)
{
  info.reset();
  double rr = NewtonResidual(r,u,f);
  bool reached = info.check(0,rr,0.);
  NewtonOutput(info);
  for(int it=1; !reached; it++)
    {
      NewtonMatrixControl(u,info);
      NewtonVectorZero(w);
      NewtonLinearSolve(w,r,info.GetLinearInfo());
      double rw = NewtonUpdate(rr,u,w,r,f,info);
      reached = info.check(it,rr,rw);
      NewtonOutput(info);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonOutput(NLInfo& nlinfo) const
{
  assert(MON!=0);
  MON->nonlinear_step(nlinfo.GetLinearInfo(),nlinfo);
}

/*-------------------------------------------------------------*/
 
void StdMultiLevelSolver::NewtonVectorZero(VectorInterface& w) const
{
  GetSolver(FinestLevel())->GetGV(w).zero();
}

/*-------------------------------------------------------------*/
 
double StdMultiLevelSolver::NewtonResidual(VectorInterface& y, const VectorInterface& x,const VectorInterface& b) const
{
  _clock_residual.start();
  DataP->CountResidual()++;
  GetSolver(ComputeLevel)->GetGV(y).equ(1.,GetSolver(ComputeLevel)->GetGV(b));
  GetSolver(ComputeLevel)->Form(y,x,-1.);
  GetSolver(ComputeLevel)->SetBoundaryVectorZero(y);
  GetSolver(ComputeLevel)->SubtractMeanAlgebraic(y);
  _clock_residual.stop();
  return NewtonNorm(y);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::NewtonMatrixControl(VectorInterface& u, NLInfo& nlinfo)
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

void StdMultiLevelSolver::NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info)
{
  _clock_solve.start();
  info.reset();
  GetSolver(FinestLevel())->GetGV(x).zero();

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
  GetSolver(ComputeLevel)->SubtractMean(x);
  _clock_solve.stop();
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Gmres(VectorInterface& x, const VectorInterface& f, CGInfo& info)
{
  int n = DataP->GmresMemSize();

  StdSolver* S = dynamic_cast<StdSolver*>(GetSolver(ComputeLevel));
  GMRES<StdSolver,StdMultiLevelSolver,VectorInterface> gmres(*S,*this,n);

  gmres.solve(x,f,info);
}

/*-------------------------------------------------------------*/

// void StdMultiLevelSolver::NewtonUpdateShowCompResiduals(VectorInterface& x, VectorInterface& r, const VectorInterface& f){

//   if( DataP->ShowCompResidualNames() ) {
//     // nur beim ersten durchgang? :  cginfo.statistics().totaliter()
//     const ComponentInformation*  CI = GetProblemDescriptor()->GetComponentInformation();
//     cout << "                                  ";
//     int     ncomps = GetSolver(FinestLevel())->GetGV(r).ncomp();
//     string  str;  
//     for(int i=0;i<ncomps;i++){
//       CI->GetScalarName(i,str);
//       printf("%-9.9s", str.c_str());
//       if(i<ncomps-1) cout << " "; 
//     }
//     cout << "  " << endl;
//   }

//   if( DataP->ShowNonLinearCompResiduals() )
//     {
//       int ncomps = GetSolver(FinestLevel())->GetGV(r).ncomp();
//       cout << "      nonlin L8-comp-residuals: [";
//       for(int i=0;i<ncomps;i++){
//         double res_comp = GetSolver(FinestLevel())->GetGV(r).CompNormL8(i)+1e-99;
//         printf(" %3.2e", res_comp);
//         if(i<ncomps-1) cout << ","; 
//       }
//       cout << " ]" << endl;
//     }

//   if( DataP->ShowLinearCompResiduals() ) 
//     { 
//       // da das residuum *moeglicherweise* nachher auch benutzt wird, mache hier einen backup 
//       GlobalVector  r_backup; 
//       GlobalVector &r_vector = r.Vector(ComputeLevel);
//       int ncomps = r_vector.ncomp(); 

//       r_backup.ncomp()  = ncomps;
//       r_backup.resize(r_vector.n());
//       r_backup.equ(1.,r_vector);
   
//       GetSolver(ComputeLevel)->MatrixResidual(r, x, f);  

//       cout << "         lin L8-comp-residuals: ["; 
//       for(int i=0;i<ncomps;i++){ 
//         double res_comp = r_vector.CompNormL8(i)+1e-99; 
//         printf(" %3.2e", res_comp); 
//         if(i<ncomps-1) cout << ",";
//       } 
//       cout << " ]" << endl; 

//       r_vector.equ(1.,r_backup);
//     } 

// }

double StdMultiLevelSolver::NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo)
{
  const CGInfo& linfo = nlinfo.GetLinearInfo();
  bool lex  = linfo.control().status()=="exploded";

  double nn = NewtonNorm(dx);
  double nr = GetSolver(ComputeLevel)->GetGV(r).norm();

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
          GetSolver(ComputeLevel)->GetGV(x).add(relax*(omega-1.),GetSolver(ComputeLevel)->GetGV(dx));
          relax *= omega;
        }
      else
        {
          GetSolver(ComputeLevel)->GetGV(x).add(relax,GetSolver(ComputeLevel)->GetGV(dx));
        }

      NewtonResidual(r,x,f);
      rr = NewtonNorm(r);
      message = nlinfo.check_damping(iter,rr);

      if (message=="ok")       break;
      if (message=="continue") continue;
      if (message=="exploded") 
        {
          GetSolver(ComputeLevel)->GetGV(x).add(-relax,GetSolver(ComputeLevel)->GetGV(dx));
          relax = 0.;
          cout << "Damping exploded !!!!!" << endl;
          nlinfo.control().status() = "diverged";
          break;
        }
    }

  //NewtonUpdateShowCompResiduals(x, r, f);

  return NewtonNorm(dx);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleMatrix(VectorInterface& u)
{
  SolutionTransfer(u);
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->MatrixZero();
      GetSolver(l)->AssembleMatrix(u,1.);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleMatrix(VectorInterface& u, NLInfo& nlinfo)
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

void StdMultiLevelSolver::ComputeIlu(VectorInterface& u)
{
  SolutionTransfer(u);
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->ComputeIlu(u);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::BoundaryInit(VectorInterface& u) const
{
  for(int l=0;l<_SP.size();l++)
    GetSolver(l)->BoundaryInit(u);
}

/*-------------------------------------------------------------*/

string StdMultiLevelSolver::LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info)
{
  ComputeLevel = level;

  GetSolver(ComputeLevel)->HNAverage(u);
  
  info.reset();
  
  int clevel=Gascoigne::max_int(DataP->CoarseLevel() ,0);
  if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

  LinearMg(ComputeLevel,clevel,u,b,info);

  GetSolver(ComputeLevel)->SubtractMean(u);

  GetSolver(ComputeLevel)->HNZero(u);
  string status = info.control().status();

  return status;
}

/*-------------------------------------------------------------*/

string StdMultiLevelSolver::Solve(int level, VectorInterface& u, const VectorInterface& b, NLInfo& nlinfo)
{
  ComputeLevel = level;

  string status;
  if(DataP->NonLinearSolve() == "newton")
    {
      GetSolver(ComputeLevel)->HNAverage(u);
      newton(u,b,_res,_cor,nlinfo);
      GetSolver(ComputeLevel)->HNZero(u);
      return nlinfo.CheckMatrix();
    }
  else
    {
      assert(0);
    }
}

/*-------------------------------------------------------------*/

double StdMultiLevelSolver::ComputeFunctional(VectorInterface& f, const VectorInterface& u, const Functional* FP) const
{
  return GetSolver(ComputeLevel)->ComputeFunctional(f,u,FP);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Transfer(int high, int low, VectorInterface& u) const
{
  for(int l=high;l>=low;l--)
    {
      Transfer(l,GetSolver(l-1)->GetGV(u),GetSolver(l)->GetGV(u));
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(int high, int low, VectorInterface& u) const
{
  for(int l=high;l>=low;l--)
    {
      SolutionTransfer(l,GetSolver(l-1)->GetGV(u),GetSolver(l)->GetGV(u));
      GetSolver(l-1)->SetBoundaryVector(u);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(VectorInterface& u) const
{
  SolutionTransfer(ComputeLevel,1,u);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const
{
  GlobalVector& uuf = const_cast<GlobalVector&>(uf);

  GetSolver(l)->GetDiscretization()->HNAverage(uuf);
  Transfer(l,ul,uf);
  GetSolver(l)->GetDiscretization()->HNZero(uuf);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const
{
  assert(_Interpolator[l-1]);
  _Interpolator[l-1]->SolutionTransfer(ul,uf);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AssembleDualMatrix(VectorInterface& u)
{
  for(int l=0; l<nlevels(); l++)
    {
      GetSolver(l)->AssembleDualMatrix(u,1.);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::DeleteVector(VectorInterface& v)
{
  for(int l=0; l<nlevels(); ++l)  
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->DeleteVector(&v);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::precondition(VectorInterface& x, VectorInterface& y)
{
  CGInfo& precinfo = DataP->GetPrecInfo();
  precinfo.reset();
  precinfo.check(0.,0.);

  int clevel=Gascoigne::max_int(DataP->CoarseLevel(),0);
  if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

  LinearMg(ComputeLevel,clevel,x,y,precinfo);
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Equ(VectorInterface& dst, double s, const VectorInterface& src) const
{
  for(int l=0; l<nlevels(); l++)
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->Equ(dst,s,src);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::Zero(VectorInterface& dst) const
{
  for(int l=0; l<nlevels(); l++)
    {
      const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(l));
      S->Zero(dst);
    }
}

/*-------------------------------------------------------------*/

void StdMultiLevelSolver::AddNodeVector(const string& name, VectorInterface& gq)
{
  Transfer(ComputeLevel,1,gq);
  for(int l=0; l<nlevels(); l++)
    {
      GetSolver(l)->AddNodeVector(name,gq);
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
