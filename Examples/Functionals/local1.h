#ifndef  __LOCAL1_h
#define  __LOCAL1_h

#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidebyequation.h"
#include  "dirichletdatabycolor.h"
#include  "pointfunctional.h"
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
    GetRightHandSidePointer() = new Gascoigne::RightHandSideByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class LocalPointFunctional : public Gascoigne::PointFunctional
{
 public:
  double J(const std::vector<double>& u) const
    {
      assert(u.size()==2);
      double a = u[1] - 0.0625;
      double b = u[0] - 0.25*0.25*0.75*0.75;
      return 0.5*(a*a + b*b);
    }

  std::string GetName() const {
    return "LeastSquaresFunctional";
  }
};

/*---------------------------------------------------*/

class LocalDragFunctional : public virtual Gascoigne::ResidualFunctional
{
 public:
  LocalDragFunctional() : ResidualFunctional()
    {
      __comps.push_back(0);
      __cols.insert(1);
      __scales.push_back(1);
      ExactValue() = 1./6.;

      __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
    }

  std::string GetName() const {
    return "LocalDrag";
  }
};

/*---------------------------------------------------*/

class LocalDomainFunctional : public virtual Gascoigne::AllDomainFunctional
{
 public:
  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 1./36.;
    }
  
  std::string GetName() const {
    return "LocalDomain";
  }
};

/*---------------------------------------------------*/

#endif

