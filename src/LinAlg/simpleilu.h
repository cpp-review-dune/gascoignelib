#ifndef  __SimpleIlu_h
#define  __SimpleIlu_h



/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleIlu

///
///
/////////////////////////////////////////////


#include  "simplematrix.h"
#include  "matrixinterface.h"

class SimpleIlu : public SimpleMatrix
{
protected:

  nvector<int>              p,q;
  mutable nvector<double>   yp;

  void hin(const nvector<double>& y) const;
  void her(nvector<double>& y) const;
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
    void solve(nvector<double>& x) const;
    void solve_transpose(nvector<double>& x) const;
};


#endif
