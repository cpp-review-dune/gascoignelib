#ifndef  __DomainMeanFunctional_h
#define  __DomainMeanFunctional_h

#include  "domainfunctional.h"
#include  <set>

/*-----------------------------------------*/

class DomainMeanFunctional  : public virtual DomainFunctional
{
protected:

  std::string  _domain;
  int          _comp;
  double       _x0, _x1, _y0, _y1, _z0, _z1;

public:

  DomainMeanFunctional(const Equation& EQ, const std::vector<std::string>& args);

  std::string GetName() const {return "domain_mean";}

  int    GetComp() const {return _comp;}

  double J(const FemFunction& U, const Vertex2d& v) const;
  double J(const FemFunction& U, const Vertex3d& v) const;
};


#endif
