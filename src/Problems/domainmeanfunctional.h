#ifndef  __DomainMeanFunctional_h
#define  __DomainMeanFunctional_h

#include  "domainfunctional.h"
#include  <set>

/*-----------------------------------------*/

class DomainMeanFunctional  : public virtual DomainFunctional
{
protected:

  std::string  _domain;
  int          _comp, _ncomp;
  double       _x0, _x1, _y0, _y1, _z0, _z1;

public:

  DomainMeanFunctional(const std::vector<std::string>& args);
  ~DomainMeanFunctional() {}

  std::string GetName() const {return "domain_mean";}

  int    GetNcomp() const {return _ncomp;}
  int    GetComp()  const {return _comp;}

  double J(const Gascoigne::FemFunction& U, const Vertex2d& v) const;
  double J(const Gascoigne::FemFunction& U, const Vertex3d& v) const;
};


#endif
