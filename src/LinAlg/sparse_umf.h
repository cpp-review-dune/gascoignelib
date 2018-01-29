/*----------------------------   sparse_umf.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __sparse_umf_H
#define __sparse_umf_H
/*----------------------------   sparse_umf.h     ---------------------------*/


#include  "iluinterface.h"
#include  "sparseblockmatrix.h"


namespace Gascoigne
{

  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments UmfIlu
  ///
  ///
  /////////////////////////////////////////////

  template<class B>
    class SparseUmf : virtual public IluInterface
  {
  protected:
    const SparseBlockMatrix<B>* __AS;


  protected:
    
    // fuer umfpack
    double *Control;
    double *Info;
    void *Symbolic, *Numeric ;

    vector<long>     __Ap,__Ac;
    nvector<double> __Ax;

    int __ncomp;
    
    
  public:

    //
    ///  Constructor 
    //
    SparseUmf(const MatrixInterface* A);
    SparseUmf();


    void SetMatrix(const MatrixInterface* A);
    
    ~SparseUmf();

    
    std::string GetName() const { return "Sparse-UMF"; }

    nvector<double>& GetRaw() { return __Ax; }
    vector<long>&    GetRawColumn() { return __Ac; }
    vector<long>&    GetRawPosition() { return __Ap; }

    
    
    
    int   n()          const { std::cerr << "SparseUmf::n()" << std::endl; abort(); }

    void zero()
    {
      __Ax.zero();
    }
    
    void ReInit(const SparseStructureInterface* SS)
    {
      assert(dynamic_cast<const SparseBlockMatrix<B> *>(__AS));
    }
    
    void copy_entries(const MatrixInterface*  A);
    void add_entries(double s, const MatrixInterface*  A);

    void ConstructStructure(const IntVector& perm, const MatrixInterface& A);
    void ConstructStructure(int ncomp, const SparseStructure& SS);
    void modify(int c, double s) {}
    void compute_ilu ();
    
    void solve(GlobalVector& x) const;
    
    void SolveTranspose(DoubleVector& x, const DoubleVector& b)
    {
      abort();
    }
    
  };

}


/*----------------------------   sparse_umf.h     ---------------------------*/
/* end of #ifndef __sparse_umf_H */
#endif
/*----------------------------   sparse_umf.h     ---------------------------*/
