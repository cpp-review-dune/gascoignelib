#include  "problemdescriptor1.h"
#include  "stdloop.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "domainmeanfunctional.h"

/*---------------------------------------------------*/

class LocalDragFunctional : public virtual ResidualFunctional
{
 public:

  LocalDragFunctional() : ResidualFunctional()
    {
      _comp = 0;
      _col.insert(1);
      _scale = 1;
      ExactValue() = 1./8.;
      beautifulname = "LocalDrag";

      _DD  = new DirichletDataByColor(GetComp(),GetColors(),GetScale());
    }
  ~LocalDragFunctional() {}
};

/*---------------------------------------------------*/

class LocalDomainFunctional : public virtual AllDomainFunctional
{
 public:

  LocalDomainFunctional() : AllDomainFunctional(1,0)
    {
      ExactValue() = 0.02776989201546093;
      beautifulname = "LocalDomain";
    }
  ~LocalDomainFunctional() {}
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Gascoigne::ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor1 LPD;
  LPD.BasicInit(&paramfile);

  /////////////
  // Loop
  /////////////
  StdLoop loop;
  loop.BasicInit(&paramfile);

  /////////////
  // Functionals
  /////////////
  LocalDragFunctional   j0; 
  LocalDomainFunctional j1;
  std::vector<const Functional*> J(2);
  J[0] = &j0;
  J[1] = &j1;
  loop.SetFunctionals(J);
  
  loop.run(&LPD);

  return 0;
}
