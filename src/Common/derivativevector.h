#ifndef  __derivativevector_h
#define  __derivativevector_h

#include  "numfixarray.h"

/*--------------------------------------------------------*/

namespace Gascoigne
{
class DerivativeVector : public numfixarray<6,double>
{
 private:
  std::map<std::string,double> _M;

 public:

  double m() const { return (*this)[0]; }
  double x() const { return (*this)[1]; }
  double y() const { return (*this)[2]; }
  double z() const { return (*this)[3]; }
  double n() const { return (*this)[4]; }
  double D() const { return (*this)[5]; }  // fuer Laplace

  double& m() { return (*this)[0]; }
  double& x() { return (*this)[1]; }
  double& y() { return (*this)[2]; }
  double& z() { return (*this)[3]; }
  double& n() { return (*this)[4]; }
  double& D() { return (*this)[5]; }

  double aux(const std::string &name) const
    {
      std::map<std::string,double>::const_iterator p = _M.find(name);
      if(p==_M.end())
        {
          std::cerr << name << " not found!" << std::endl;
          abort();
        }
      return p->second;
    }

  double &aux(const std::string &name) { return _M[name]; }

};
}

/*--------------------------------------------------------*/

#endif
