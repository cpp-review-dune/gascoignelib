#ifndef  __LOCAL1_h
#define  __LOCAL1_h

#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "dirichletdatabycolor.h"
#include  "newpointfunctional.h"
#include  "residualfunctional.h"
#include  "domainmeanfunctional.h"

/*---------------------------------------------------*/

class PolynomialExactSolution : public Gascoigne::ExactSolution
{
 public:
  std::string GetName() const {return "PolynomialExactSolution";}
  int GetNcomp() const { return 1; }
  double operator()(int c, const Gascoigne::Vertex2d& v) const{
    return v.x()*(1.-v.x())*v.y()*(1.-v.y());
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
    GetRightHandSideDataPointer() = new Gascoigne::RightHandSideDataByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class LocalPointFunctional : public Gascoigne::NewPointFunctional
{
 public:
 LocalPointFunctional() : Gascoigne::NewPointFunctional() {beautifulname = "LeastSquaresFunctional";}

  double J(const std::vector<double>& u) const
    {
      assert(u.size()==2);
      double a = u[1] - 0.0625;
      double b = u[0] - 0.25*0.25*0.75*0.75;
      return 0.5*(a*a + b*b);
    }
};

/*---------------------------------------------------*/

class LocalDragFunctional : public virtual Gascoigne::ResidualFunctional
{
 public:
  LocalDragFunctional() : ResidualFunctional()
    {
      _comp = 0;
      _col.insert(1);
      _scale = 1;
      ExactValue() = 1./6.;
      beautifulname = "LocalDrag";

      _DD  = new Gascoigne::DirichletDataByColor(GetComp(),GetColors(),GetScale());
    }
};

/*---------------------------------------------------*/

class LocalDomainFunctional : public virtual Gascoigne::AllDomainFunctional
{
 public:
  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 1./36.;
      beautifulname = "LocalDomain";
    }
};

/*---------------------------------------------------*/

#endif

