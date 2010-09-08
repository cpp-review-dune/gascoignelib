#ifndef __givensrotation_h
#define __givensrotation_h

#include "solverinterface.h"

/*---------------------------------------------------------------*/
 
namespace Gascoigne
{
class GivensRotation
{
  //
  // Data
  //
  int                  n;
  DoubleMatrix  H;
  DoubleVector  ci, si,gamma;

 public:
  
  GivensRotation(int nn, double) ;
  double&         matrix(int i, int j)             { return H(i,j);}
  void               givens(double& a, double& b, int i) const;
  double           orthogonalization(int dim) ;
  DoubleVector getcoefficients();
};
} 

/*---------------------------------------------------------------*/
 
#endif
