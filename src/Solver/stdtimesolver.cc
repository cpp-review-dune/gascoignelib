#include  "stdtimesolver.h"
#include  "simplematrix.h"
#include  "compose_name.h"

using namespace std;

/*-------------------------------------------------------------*/
  
StdTimeSolver::StdTimeSolver()
  : StdSolver(), _MMP(NULL), dt(0.), theta(0.), time(0.)
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
  
void StdTimeSolver::SetTimeData(double d, double th, double ti) 
{
  dt    = d;
  theta = th;
  time  = ti;

  GetProblemDescriptor()->SetTime(time,dt);
}

/*-------------------------------------------------------------*/

void StdTimeSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  const Equation* EQ = PDX.GetEquation();
  assert(EQ);
  EQ->SetTimePattern(GetTimePattern());

  int ncomp = EQ->ncomp();
  if (GetMassMatrixPointer()==NULL)
    GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
  
  StdSolver::SetProblem(PDX);
}

/*-------------------------------------------------------------*/

MatrixInterface* StdTimeSolver::NewMassMatrix(int ncomp, const string& matrixtype)
{
  return new SimpleMatrix;
}

/*-------------------------------------------------------*/

void StdTimeSolver::MemoryMatrix()
{
  ConstructPressureFilter();

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
  
void StdTimeSolver::BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)
{
  StdSolver::BasicInit(level, paramfile, MP);
}

/*-------------------------------------------------------*/

void StdTimeSolver::TimeRhs(BasicGhostVector& gf, const BasicGhostVector& gu) const
{
  assert(theta>0.);

  double d = -(1.-theta)/theta;

  StdSolver::Residual(gf,gu,d);

  GlobalVector& f = GetGV(gf);
  const GlobalVector& u = GetGV(gu);

  d = 1./(dt*theta);

  GetMassMatrix()->vmult_time(f,u,GetTimePattern(),d);

  StdSolver::Rhs(gf,1./theta);
}

/*-------------------------------------------------------*/

void StdTimeSolver::Residual(BasicGhostVector& gy, const BasicGhostVector& gx, double d) const
{
  StdSolver::Residual(gy,gx,d);

  if (dt==0.) return;
  assert(theta>0.);

  double scale = d/(dt*theta);

  const GlobalVector& x = GetGV(gx);
  GlobalVector& y = GetGV(gy);
  GetMassMatrix()->vmult_time(y,x,GetTimePattern(),scale);
}

/*-------------------------------------------------------*/

void StdTimeSolver::AssembleMatrix(BasicGhostVector& gu, double d)
{
  StdSolver::AssembleMatrix(gu,d);

  if (dt==0.) return;
  assert(theta>0.);

  double scale = 1./(dt*theta);
  GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),GetTimePattern(),scale);

  StdSolver::DirichletMatrix();
}

/*-------------------------------------------------------*/

void StdTimeSolver::RhsL2Projection(BasicGhostVector& f, const BasicGhostVector& u) const
{
  HNAverage(u); 
  const InitialCondition* IC  = GetProblemDescriptor()->GetInitialCondition();
  if(IC==NULL) return;

  GetMeshInterpretor()->Rhs(GetGV(f),*IC,1.);

  HNDistribute(f);
  HNZero(u);
}

/*-------------------------------------------------------*/

void StdTimeSolver::L2Projection(BasicGhostVector& Gu, const BasicGhostVector& Gf)
{
  BasicGhostVector Gg, Gr, Gd;

  Gg.SetName("g");
  Gr.SetName("r");
  Gd.SetName("d");
  RegisterVector(Gg);
  RegisterVector(Gr);
  RegisterVector(Gd);

  MemoryVector();

  GlobalVector& u = GetGV(Gu);
  const GlobalVector& f = GetGV(Gf);
  GlobalVector& g = GetGV(Gg);
  GlobalVector& r = GetGV(Gr);
  GlobalVector& d = GetGV(Gd);

  int    MaxIter = 100;
  double Tol     = 1e-8;

  u.zero();
  if (f.norm()<1.e-16) return;
  // ??????????? geandert: war vorher
//   r.equ(-1,f);
//   d.equ(-1,f);
  // was zur loesung -u fuehrte
  r.equ(1,f);
  d.equ(1,f);
  double Res = r*r;
  double FirstRes = Res;
  cerr << "\t\tcg " << 0 << "\t" << sqrt(FirstRes) << endl;

  for(int iter=1;iter<=MaxIter;iter++)
    {
      g.zero();
      GetMassMatrix()->vmult_time(g,d,GetTimePattern());
      double lambda = Res/(g*d);

      u.add(lambda,d);
      r.add(-lambda,g);

      double ResOld = Res;
      Res = r*r;
      cerr << "\t\tcg " << iter << "\t" << sqrt(Res) << "\n";
      if (Res<Tol*Tol*FirstRes) 
	{
	  return;
	}
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }
}
