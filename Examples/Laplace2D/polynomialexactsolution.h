#ifndef  __PolynomialExactSolution_h
#define  __PolynomialExactSolution_h

#include  "exactsolution.h"

/////////////////////////////////////////////
///
///@brief
///  ... comments PolynomialExactSolution

///
///
/////////////////////////////////////////////




class PolynomialExactSolution : public ExactSolution
{
public:


//
///  Constructor 
//
  PolynomialExactSolution() : ExactSolution() {}

  std::string GetName() const {return "PolynomialExactSolution";}
  double operator()(int c, const Vertex2d& v)const ;
  int GetNcomp() const { return 1; }

};


#endif
