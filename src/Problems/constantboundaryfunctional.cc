#include  "constantboundaryfunctional.h"


/*-----------------------------------------*/

ConstantBoundaryFunctional::ConstantBoundaryFunctional()
{
}

/*-----------------------------------------*/

ConstantBoundaryFunctional::~ConstantBoundaryFunctional()
{
}

/*-----------------------------------------*/

ConstantBoundaryFunctional::ConstantBoundaryFunctional
(const std::vector<std::string>& args)
{
  Construct(args);
}

/*-----------------------------------------*/

void ConstantBoundaryFunctional::Construct
(const std::vector<std::string>& args)
{
  // Baustelle
  AddColor(atoi(args[0].c_str()));
  comp  = atoi(args[1].c_str());
  value = atof(args[2].c_str());
}

/*-----------------------------------------*/

double ConstantBoundaryFunctional::J
(const FemFunction& U, const Vertex2d& v) const
{
  return value*U[comp].m();
}
