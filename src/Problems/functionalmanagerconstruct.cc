#include  "functionalmanager.h"
#include  "boundaryfunctional.h"
#include  "domainmeanfunctional.h"
#include  "residualfunctional.h"
#include  "pointfunctional.h"
#include  "zerofunctional.h"
#include  "constantboundaryfunctional.h"
#include  "stringutil.h"

/*-----------------------------------------*/

Functional* FunctionalManager::ConstructFunctional
(const std::string& name, const std::string& param, const Equation& EQ)
{
  std::vector<std::string> args = StringSplit(param.c_str(),'_');

  if(name=="zero")
    {
      return new ZeroFunctional;
    }
  if(name=="residual")
    {
      return new ResidualFunctional(EQ, args);
    }
  if(name=="point")
    {
      return new PointFunctional(EQ,args);
    }
  if(name=="domainmean")
    {
      return new DomainMeanFunctional(EQ,args);
    }
  if(name=="constantboundary")
    {
      return new ConstantBoundaryFunctional(EQ,args);
    }
  if(name=="nusselt")
    {
      assert(0);
    }
  std::cerr << "FunctionalManager::ConstructFunctional()\n";
  std::cerr << "unknown functional name: " << name << std::endl;
  abort();
}
