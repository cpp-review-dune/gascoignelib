#ifndef  __DomainMeanFunctional_h
#define  __DomainMeanFunctional_h

#include  "domainfunctional.h"
#include  <set>

/*-----------------------------------------*/

namespace Gascoigne
{
class AllDomainFunctional  : public virtual DomainFunctional
{
protected:

  int  _comp, _ncomp;

public:

  AllDomainFunctional(int nc, int c) { _ncomp = nc; _comp = c; }
  ~AllDomainFunctional() {}

  std::string GetName() const {return "AllDomainFunctional";}

  int    GetNcomp() const {return _ncomp;}
  int    GetComp()  const {return _comp;}

  double J(const FemFunction& U, const Vertex2d& v) const;
  double J(const FemFunction& U, const Vertex3d& v) const;
};

/*-----------------------------------------*/

class SubDomainFunctional  : public AllDomainFunctional
{
protected:

  double  _x0, _x1, _y0, _y1, _z0, _z1;

public:

  SubDomainFunctional(int nc, int c) : AllDomainFunctional(nc,c) {};
  ~SubDomainFunctional() {}

  std::string GetName() const {return "SubDomainFunctional";}

  void SetCoordinates(double x0, double x1, double y0, double y1);

  double J(const FemFunction& U, const Vertex2d& v) const;
  double J(const FemFunction& U, const Vertex3d& v) const;
};
}

#endif
