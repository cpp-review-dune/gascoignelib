#ifndef  __Laplace3d_h
#define  __Laplace3d_h

#include  "laplace2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class Laplace3d : public Laplace2d
{
  double betax, betay, betaz;

public:

  Laplace3d(const ParamFile* pf);

  std::string GetName() const { return "Laplace3d";}

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& D, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  void OperatorStrong(DoubleVector& b, const FemFunction& U)const;
};
}

#endif
