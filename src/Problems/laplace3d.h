#ifndef  __Laplace3d_h
#define  __Laplace3d_h

#include  "laplace2d.h"

/*-----------------------------------------*/

class Laplace3d : public Laplace2d
{
  double betax, betay, betaz;

public:

  Laplace3d(const Gascoigne::ParamFile* pf);

  std::string GetName() const { return "Laplace3d";}

  //
  // Semilinear Form
  //

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& D, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};


#endif
