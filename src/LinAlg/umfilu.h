#ifndef  __UmfIlu_h
#define  __UmfIlu_h

#ifdef __WITH_UMFPACK__

#include  "iluinterface.h"
#include  "simplematrix.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments UmfIlu
///
///
/////////////////////////////////////////////

class UmfIlu : virtual public IluInterface, public SimpleMatrix
{
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
    
    std::string GetName() const { return "UmfIlu"; }
    
    void ReInit(const SparseStructureInterface* SS);

    int   n()         const { return GetStencil()->n();};

    void copy_entries(const MatrixInterface&  A);
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A);
    void Factorize();
    void Solve(DoubleVector& x, const DoubleVector& b);
    void SolveTranspose(DoubleVector& x, const DoubleVector& b);
};
}

#endif
#endif
