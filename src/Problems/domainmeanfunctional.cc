#include  "domainmeanfunctional.h"
#include  "zerodirichletdata.h"
#include  "constantrighthandside.h"


using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

double AllDomainFunctional::J(const FemFunction& U, const Vertex2d& V) const
{
  return U[_comp].m();
}

/*-----------------------------------------*/

double AllDomainFunctional::J(const FemFunction& U, const Vertex3d& V) const
{
  return U[_comp].m();
}
/*-----------------------------------------*/

void SubDomainFunctional::SetCoordinates(double x0, double x1, double y0, double y1)
{
  _x0 = x0;
  _x1 = x1;
  _y0 = y0;
  _y1 = y1;
}

/*-----------------------------------------*/

double SubDomainFunctional::J(const FemFunction& U, const Vertex2d& V) const
{
  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;

  return U[_comp].m();
}

/*-----------------------------------------*/

double SubDomainFunctional::J(const FemFunction& U, const Vertex3d& V) const
{
  if ((V.x()>_x1) || (V.x()<_x0)) return 0.;
  if ((V.y()>_y1) || (V.y()<_y0)) return 0.;
  if ((V.z()>_z1) || (V.z()<_z0)) return 0.;

  return U[_comp].m();
}
