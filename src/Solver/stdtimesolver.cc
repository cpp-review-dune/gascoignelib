#include  "stdtimesolver.h"
#include  "simplematrix.h"
#include  "compose_name.h"

using namespace std;
using namespace Gascoigne;

/*-------------------------------------------------------------*/
  
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
  
void StdTimeSolver::SetTimeData(double dt, double theta, double time, double rhs) 
{
  _dt    = dt;
  _theta = theta;
  _time  = time;
//   _rhs   = (rh<0) ? (1.-theta)/theta : rh;
  _rhs = rhs;
  if(rhs<0)  _rhs = (1.-theta)/theta;

  GetProblemDescriptor()->SetTime(_time,_dt);
}

/*-------------------------------------------------------------*/

void StdTimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();
  assert(EQ);
  EQ->SetTimePattern(GetTimePattern());
  
  StdSolver::SetProblem(PDX);
}

/*-------------------------------------------------------*/

void StdTimeSolver::RegisterMatrix()
{
  const Equation*  EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  int ncomp = EQ->ncomp();

  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
  
  StdSolver::RegisterMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::ReInitMatrix() 
{
  GetMeshInterpretor()->InitFilter(PF);
  SparseStructure SA;
  GetMeshInterpretor()->Structure(&SA);

  GetMatrix()->ReInit(&SA);
  GetIlu()->ReInit(&SA);
  GetMassMatrix()->ReInit(&SA);

  GetMassMatrix()->zero();
  GetMeshInterpretor()->MassMatrix(*GetMassMatrix());  

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

void StdTimeSolver::TimeRhsOperator(BasicGhostVector& gf, const BasicGhostVector& gu) const
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

void StdTimeSolver::TimeRhs(int k, BasicGhostVector& gf) const
{
  StdSolver::Rhs(gf,_rhs);
}

/*-------------------------------------------------------*/

void StdTimeSolver::Form(BasicGhostVector& gy, const BasicGhostVector& gx, double d) const
{
  StdSolver::Form(gy,gx,d);

  if (_dt==0.) return;
  assert(_theta>0.);

  double scale = d/(_dt*_theta);

  const GlobalVector& x = GetGV(gx);
  GlobalVector& y = GetGV(gy);
  GetMassMatrix()->vmult_time(y,x,GetTimePattern(),scale);
}

/*-------------------------------------------------------*/

void StdTimeSolver::AssembleMatrix(BasicGhostVector& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);

  if (_dt==0.) return;
  assert(_theta>0.);

  double scale = 1./(_dt*_theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),scale);

  StdSolver::DirichletMatrix();
}



/*-------------------------------------------------------*/

void StdTimeSolver::L2Projection(BasicGhostVector& Gu, const InitialCondition* IC)
{
  if(IC==NULL) return;

  GlobalVector& u = GetGV(Gu);

  BasicGhostVector Gg, Gr, Gd, Gf;

  Gg.SetName("g");
  Gr.SetName("r");
  Gd.SetName("d");
  Gd.SetName("f");
  RegisterVector(Gg);
  RegisterVector(Gr);
  RegisterVector(Gd);
  RegisterVector(Gf);

  ReInitVector();
  GlobalVector& g = GetGV(Gg);
  GlobalVector& r = GetGV(Gr);
  GlobalVector& d = GetGV(Gd);
  GlobalVector& f = GetGV(Gf);

  assert(u.ncomp()==g.ncomp());
  assert(u.n()==g.n());

  f.zero();
  GetMeshInterpretor()->Rhs(f,*IC,1.);
  HNDistribute(f);

  int    MaxIter = 100;
  double Tol     = 1e-8;

  u.zero();
  if (f.norm()<1.e-16) return;

  r.equ(1,f);
  d.equ(1,f);
  double Res = r*r;
  double FirstRes = Res;
  cout << "\t\tcg " << 0 << "\t" << sqrt(FirstRes) << endl;

  TimePattern TP(u.ncomp());
  TP.zero();
  for (int i=0; i<u.ncomp(); i++) TP(i,i) = 1.;

  for(int iter=1;iter<=MaxIter;iter++)
    {
      g.zero();
      GetMassMatrix()->vmult_time(g,d,TP);
      double lambda = Res/(g*d);

      u.add(lambda,d);
      r.add(-lambda,g);

      Res = r*r;
      cout << "\t\tcg " << iter << "\t" << sqrt(Res) << "\n";
      if (Res<Tol*Tol*FirstRes) 
	{
	  return;
	}
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }
}
