#ifndef  __LOCAL1_h
#define  __LOCAL1_h

#include  "problemdescriptorbase.h"
#include  "laplace3d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"

/*---------------------------------------------------*/

class PolynomialExactSolution : public Gascoigne::ExactSolution
{
  double quadratic(double x) const { return x*(1.-x);}
 public:
  std::string GetName() const {return "PolynomialExactSolution";}
  int GetNcomp() const { return 1; }
  double operator()(int c, const Gascoigne::Vertex3d& v) const{
    return quadratic(v.x()) * quadratic(v.y()) * quadratic(v.z());
  }
};

class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Local";}
  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new Gascoigne::Laplace3d(GetParamFile());
    GetExactSolutionPointer() = new PolynomialExactSolution();
    GetRightHandSideDataPointer() = new Gascoigne::RightHandSideDataByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

class LocalDragFunctional : public virtual Gascoigne::ResidualFunctional
{
 public:
  LocalDragFunctional() : ResidualFunctional()
    {
      _comp = 0;
      _col.insert(9);
      _scale = 1;
      ExactValue() = 1./8.;
      beautifulname = "LocalDrag";

      _DD  = new Gascoigne::DirichletDataByColor(GetComp(),GetColors(),GetScale());
    }
};

class LocalDomainFunctional : public virtual Gascoigne::AllDomainFunctional
{
 public:
  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 0.02776989201546093;
      beautifulname = "LocalDomain";
    }
};


#endif

