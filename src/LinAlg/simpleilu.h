#ifndef  __SimpleIlu_h
#define  __SimpleIlu_h

#include  "simplematrix.h"
#include  "matrixinterface.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleIlu

///
///
/////////////////////////////////////////////

class SimpleIlu : public SimpleMatrix
{
protected:

  IntVector              p,q;
  mutable DoubleVector   yp;

  void hin(const DoubleVector& y) const;
  void her(DoubleVector& y) const;
  void backward() const;
  void forward () const;
  void backward_transpose() const;
  void forward_transpose () const;

public:

//
///  Constructor 
//

    SimpleIlu() : SimpleMatrix() {}
    
    void zero() {SimpleMatrix::zero();}
    void ReInit(int n, int nentries);
    void copy_entries(const MatrixInterface*  A);
    void compute_ilu();
    void solve(DoubleVector& x) const;
    void solve_transpose(DoubleVector& x) const;
};
}

#endif
