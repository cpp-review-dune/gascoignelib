#ifndef  __LOCAL_h
#define  __LOCAL_h

#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidebyequation.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"
#include  "stdloop.h"
#include  "meshinterpolator.h"

/*---------------------------------------------------*/

class PolynomialExactSolution : public Gascoigne::ExactSolution
{
 public:
  std::string GetName() const {return "PolynomialExactSolution";}
  int GetNcomp() const { return 1; }
  double operator()(int c, const Gascoigne::Vertex2d& v) const{
    //   return v.x()*(1.-v.x())*v.y()*(1.-v.y());
    return v.x()*(1.-v.x())*v.y()*(1.-v.y());// *  (exp(v.x()+v.y()));
  }
};

// for use with slit.param !!
class SlitExactSolution : public Gascoigne::ExactSolution
{
public:
  double operator()(int c, const Gascoigne::Vertex2d& v)const 
  {
    double x = v.x();
    double y = v.y();
    double r = sqrt(x*x+y*y);
    
    double pi = Gascoigne::pi();
    double theta;

    double fx = fabs(x);
    double fy = fabs(y);
    if(fx)
      {
        theta = atan(fy/fx);

        if     ( (x<0)&&(y>=0)) theta = pi-theta;
        else if( (x<0)&&(y<0))  theta += pi;
        else if( (x>0)&&(y<0))  theta = 2.*pi-theta;
      }
    else
      {
        if(y>=0) theta = 0.5*pi;
        else     theta = 1.5*pi;
      }
    return pow(r,0.5)*sin(0.5*theta);
  }
};

/*---------------------------------------------------*/

class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Local";}
  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new Gascoigne::Laplace2d;
    GetExactSolutionPointer() = new PolynomialExactSolution();
    GetRightHandSidePointer() = new Gascoigne::RightHandSideByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class ProjectionEquation : public Gascoigne::Equation
{
 public:

  ProjectionEquation() : Gascoigne::Equation() { }
  ~ProjectionEquation() { }

  int GetNcomp() const { return 1; }
  std::string GetName() const { return "ProjectionEquation"; }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const
  {
    b[0] += U[0].m() * N.m();
  }
  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const
  {
    A(0,0) += M.m() * N.m();
  }
};

/*---------------------------------------------------*/

class ProjectionRightHandSide : public Gascoigne::DomainRightHandSide
{
 protected:
   mutable Gascoigne::FemFunction _U;
   int _ncomp;

 public:

  ProjectionRightHandSide(const Gascoigne::Equation* EQ) : Gascoigne::DomainRightHandSide() { _ncomp = EQ->GetNcomp(); }
  ~ProjectionRightHandSide() { }
  
  int GetNcomp() const { return _ncomp; }
  std::string GetName() const { return "ProjectionRightHandSide"; }

  void SetFemData(Gascoigne::FemData& q) const
  {
    assert(q.count("U")==1);
    _U = q["U"];
  }

  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Gascoigne::Vertex2d& v) const 
  {
    b[0] += _U[0].m() * N.m();
  }
};

/*---------------------------------------------------*/

class ProjectionProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Projection";}
  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new ProjectionEquation();
    GetRightHandSidePointer() = new ProjectionRightHandSide(GetEquation());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class LocalLoop : public Gascoigne::StdLoop
{
 public:

  LocalLoop() : StdLoop() { }
  void run(const Gascoigne::ProblemDescriptorInterface* PD) {
    Gascoigne::MultiLevelGhostVector u("u"), f("f");

    GetMultiLevelSolver()->RegisterVector(u);
    GetMultiLevelSolver()->RegisterVector(f);
    
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    GetMultiLevelSolver()->ReInit(*PD);
    GetMultiLevelSolver()->GetSolver()->OutputSettings();

    u.zero();
    f.zero();
    Gascoigne::MeshInterpolator MI;
    MI.BasicInit(GetMultiLevelSolver()->GetSolver(),GetMeshAgent(),"start2",PD);
    MI.RhsForProjection(f);

    GetMultiLevelSolver()->ReInit(*PD);
    GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
    GetMultiLevelSolver()->GetSolver()->Visu("ziel2",u,0);
  }
};

#endif

