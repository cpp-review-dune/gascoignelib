#include  "stdtimesolver.h"
#include  "simplematrix.h"
#include  "compose_name.h"

using namespace std;

/*-------------------------------------------------------------*/
  
namespace Gascoigne
{
StdTimeSolver::StdTimeSolver()
  : StdSolver(), _MMP(NULL), _dt(0.), _theta(0.), _time(0.), _rhs(0.)
{}

/*-------------------------------------------------------------*/
  
StdTimeSolver::~StdTimeSolver()
{
  if (_MMP) { delete _MMP; _MMP=NULL;}
}

/*-------------------------------------------------------*/

string StdTimeSolver::GetName() const
{
  return "StdTimeSolver";
}

/*-------------------------------------------------------------*/
  
void StdTimeSolver::SetTimeData(double dt, double theta, double time, double oldrhs, double newrhs) 
{
  _dt     = dt;
  _theta  = theta;
  _time   = time;
  _rhs[0] = oldrhs;
  _rhs[1] = newrhs;
  if(oldrhs<0)
  {
    _rhs[0] = (1.-theta)/theta;
    _rhs[1] = 1.;
  }

  GetProblemDescriptor()->SetTime(_time,_dt);
}

/*-------------------------------------------------------------*/

void StdTimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();

  if (EQ) 
    {
      GetTimePattern().reservesize(EQ->GetNcomp(),EQ->GetNcomp(),0.);
      EQ->SetTimePattern(GetTimePattern());
    }
  
  StdSolver::SetProblem(PDX);
}

/*-------------------------------------------------------*/

void StdTimeSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->GetNcomp();

  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
  
  StdSolver::RegisterMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::ReInitMatrix() 
{
  GetDiscretization()->InitFilter(_PF);
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);

  GetMassMatrix()->ReInit(&SA);
  GetMassMatrix()->zero();
  GetDiscretization()->MassMatrix(*GetMassMatrix());  

//   string name("masse");
//   compose_name(name,mylevel);
//   ofstream file(name.c_str());
//   GetMassMatrix()->Write(file);
}

/*-------------------------------------------------------------*/

MatrixInterface* StdTimeSolver::NewMassMatrix(int ncomp, const string& matrixtype)
{
  return new SimpleMatrix;
}

/*-------------------------------------------------------------*/
  
void StdTimeSolver::BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)
{
  StdSolver::BasicInit(level, paramfile, MP);
}

/*-------------------------------------------------------*/

void Gascoigne::StdTimeSolver::InitialCondition(VectorInterface& f, double d) const
{
  StdTimeSolver::InitialCondition(GetGV(f),d);
}

/*-------------------------------------------------------*/

void Gascoigne::StdTimeSolver::InitialCondition(GlobalVector& f, double d) const
{
  HNAverageData();

  const Application* IC  = GetProblemDescriptor()->GetInitialCondition();

  if(IC)
    {
       bool done=false;
       const DomainInitialCondition *DRHS = dynamic_cast<const DomainRightHandSide *>(IC);
       if(DRHS)
       {
         GetDiscretization()->Rhs(f,*DRHS,d);
         done = true;
       }
       const DiracInitialCondition *NDRHS = dynamic_cast<const DiracRightHandSide *>(IC);
       if(NDRHS)
       {
         GetDiscretization()->DiracRhs(f,*NDRHS,d);
         done =true;
       }
       if(!done)
       {
         cerr << "InitialCondition should be either of type DomainRightHandSide or DiracRightHandSide!!!" << endl;
         abort();
       }
    }

  HNZeroData();
  if (distribute())
  {
    GetDiscretization()->HNDistribute(f);
  }
}

/*-------------------------------------------------------*/

void StdTimeSolver::TimeRhsOperator(VectorInterface& gf, const VectorInterface& gu) const
{
  assert(_theta>0.);
  double d = -(1.-_theta)/_theta;
  StdSolver::Form(gf,gu,d);

  if (_dt>0.)
    {
      GlobalVector& f = GetGV(gf);
      const GlobalVector& u = GetGV(gu);
      
      double d = 1./(_dt*_theta);
      GetMassMatrix()->vmult_time(f,u,GetTimePattern(),d);
    }
}

/*-------------------------------------------------------*/

void StdTimeSolver::TimeRhs(int k, VectorInterface& gf) const
{
  StdSolver::Rhs(gf,_rhs[k-1]);
}

/*-------------------------------------------------------*/

void StdTimeSolver::Form(VectorInterface& gy, const VectorInterface& gx, double d) const
{
  StdSolver::Form(gy,gx,d);

  if (_dt==0.) return;
  assert(_theta>0.);

  double scale = d/(_dt*_theta);

  const GlobalVector& x = GetGV(gx);
  GlobalVector& y = GetGV(gy);

  assert(y.n()==x.n());
  GetMassMatrix()->vmult_time(y,x,GetTimePattern(),scale);
}

/*-------------------------------------------------------*/

void StdTimeSolver::AssembleMatrix(const VectorInterface& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);

  if (_dt==0.) return;
  assert(_theta>0.);

  double scale = d/(_dt*_theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),scale);

  StdSolver::DirichletMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::L2Projection(VectorInterface& Gu)
{
  GlobalVector& u = GetGV(Gu);

  GlobalVector f(u.ncomp(),u.n());

  TimePattern TP(u.ncomp());
  TP.zero();
  for (int i=0; i<u.ncomp(); i++) TP(i,i) = 1.;

  u.zero();
  f.zero();

  InitialCondition(f);

  PrecondCGMass(u,f,TP);
}

/*-------------------------------------------------------*/

string StdTimeSolver::PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s)
{
  bool reached;
  int iter = 0;
  
  GlobalVector g(u.ncomp(),u.n());
  GlobalVector r(u.ncomp(),u.n());
  GlobalVector d(u.ncomp(),u.n());

  assert(u.ncomp()==g.ncomp());
  assert(u.n()==g.n());

  SimpleMatrix *SM = dynamic_cast<SimpleMatrix *>(GetMassMatrix());
  assert(SM);
  SM->PrepareJacobi(s);

  SM->vmult_time(f,u,TP,-s);

  SM->JacobiVector(u);
  SM->JacobiVectorInv(f);

  r.equ(1,f);
  d.equ(1,f);
  double Res = r*r;
  double FirstRes = Res;
  cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
  
  if (sqrt(Res)<_Dat.GetCgMassGlobalTol()) 
  {
    reached = true;
  }
  else
  {
    reached = false;
  }

  while(!reached && iter<_Dat.GetCgMassMaxIter())
    {
      iter++;
      g.zero();
      SM->vmult_time_Jacobi(g,d,TP,s);
      double lambda = Res/(g*d);

      u.add(lambda,d);
      r.add(-lambda,g);

      Res = r*r;
      cout << "\t\tpcg " << iter << "\t" << sqrt(Res) << endl;
      if (Res < _Dat.GetCgMassTol() * _Dat.GetCgMassTol() * FirstRes || sqrt(Res)<_Dat.GetCgMassGlobalTol()) 
      {
        reached = true;
      }
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }

  SM->JacobiVectorInv(u);

  if(iter==_Dat.GetCgMassMaxIter())
  {
    return "too many iterations";
  }
  else
  {
    return "converged";
  }
}

/*-------------------------------------------------------*/

void StdTimeSolver::SetMassMatrix(Gascoigne::MatrixInterface &MM, bool init)
{
  if(init)
  {
    SparseStructure SA;
    GetDiscretization()->Structure(&SA);

    MM.ReInit(&SA);
    MM.zero();
    GetDiscretization()->MassMatrix(MM);
  }
  GetMassMatrixPointer() = &MM;
}
}
