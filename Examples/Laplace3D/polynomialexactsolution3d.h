
#ifndef  __PolynomialExactSolution3d_h
#define  __PolynomialExactSolution3d_h

#include  "exactsolution.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments PolynomialExactSolution
///
///
/////////////////////////////////////////////


class PolynomialExactSolution3d : public ExactSolution
{
  double quadratic(double x) const { return x*(1.-x);}

public:

//
///  Constructor 
//
  PolynomialExactSolution3d() : ExactSolution() {}

  std::string GetName() const {return "PolynomialExactSolution3d";}
  double operator()(int c, const Vertex3d& v)const ;
  int GetNcomp() const { return 1; }
};


#endif
