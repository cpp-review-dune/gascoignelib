#ifndef  __UmfIlu_h
#define  __UmfIlu_h



/////////////////////////////////////////////
///
///@brief
///  ... comments UmfIlu

///
///
/////////////////////////////////////////////


#include  "iluinterface.h"
#include  "simplematrix.h"


class UmfIlu : virtual public IluInterface, public SimpleMatrix
{
public:


private:

  const SimpleMatrix* AP;

protected:

  // fuer umfpack
  double *Control;
  double *Info;
  void *Symbolic, *Numeric ;

public:


  //
  ///  Constructor 
    //
    UmfIlu(const MatrixInterface* A);
    ~UmfIlu();
    
    void ReInit(const SparseStructureInterface* SS);
    void copy_entries(const MatrixInterface&  A);
    void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A);
    void Factorize();
    void Solve(nvector<double>& x, const nvector<double>& b);
    void SolveTranspose(nvector<double>& x, const nvector<double>& b);
};


#endif
