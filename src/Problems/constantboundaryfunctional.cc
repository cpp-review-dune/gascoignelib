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
(const Equation& EQ, const std::vector<std::string>& args)
{
  Construct(EQ,args);
}

/*-----------------------------------------*/

void ConstantBoundaryFunctional::Construct
(const Equation& EQ, const std::vector<std::string>& args)
{
  int ncomp = EQ.ncomp();
  
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
