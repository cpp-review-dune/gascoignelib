#ifndef  __derivativevector_h
#define  __derivativevector_h

#include  "numfixarray.h"

/*--------------------------------------------------------*/

class DerivativeVector : public numfixarray<6,double>
{
 public:

  ~DerivativeVector() {};

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
};

/*--------------------------------------------------------*/

#endif
