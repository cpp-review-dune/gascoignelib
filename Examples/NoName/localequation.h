#ifndef  __LocalEquation_h
#define  __LocalEquation_h

#include  "navierstokesgls2d.h"

/////////////////////////////////////////////
////
////@brief
////  ... comments LocalEquation

////
////
/////////////////////////////////////////////



class LocalEquation : public NavierStokesGls2d
{
private:


protected:

public:


//
////  Con(De)structor 
//

  LocalEquation(const Gascoigne::ParamFile* paramfile) : NavierStokesGls2d(paramfile) {}
  ~LocalEquation() {}

  void point(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const;

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;

};


#endif
