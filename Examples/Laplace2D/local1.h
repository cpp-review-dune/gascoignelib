#ifndef  __LOCAL1_h
#define  __LOCAL1_h

#include  "problemdescriptorbase.h"
#include  "laplace2d.h"
#include  "dirichletdatabyexactsolution.h"
#include  "righthandsidedatabyequation.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"

/*---------------------------------------------------*/

class PolynomialExactSolution : public Gascoigne::ExactSolution
{
 public:
  std::string GetName() const {return "PolynomialExactSolution";}
  int GetNcomp() const { return 1; }
  double operator()(int c, const Gascoigne::Vertex2d& v) const{
    //   return v.x()*(1.-v.x())*v.y()*(1.-v.y());
    return v.x()*(1.-v.x())*v.y()*(1.-v.y()) *  (exp(v.x()+v.y()));
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

class ProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Local";}
  void BasicInit(const Gascoigne::ParamFile* pf) {
    GetEquationPointer() = new Gascoigne::Laplace2d;
    GetExactSolutionPointer() = new PolynomialExactSolution();
    GetRightHandSideDataPointer() = new Gascoigne::RightHandSideDataByEquation(GetEquation(), GetExactSolution());
    GetDirichletDataPointer() = new Gascoigne::DirichletDataByExactSolution(GetExactSolution());
    GetBoundaryManagerPointer() = new Gascoigne::BoundaryManager(pf);
    
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
      _col.insert(1);
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

