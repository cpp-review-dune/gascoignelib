#ifndef  __Laplace3d_h
#define  __Laplace3d_h

#include  "laplace.h"

/*-----------------------------------------*/

class Laplace3d : public Laplace
{
  double betax, betay, betaz;

public:

  Laplace3d(const ParamFile* pf);

  std::string GetName() const { return "Laplace3d";}

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& D, const FemFunction& U, const DerivativeVector& M, const TestFunction& N) const;
};


#endif
