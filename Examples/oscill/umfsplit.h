/*----------------------------   umfsplit.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __umfsplit_H
#define __umfsplit_H
/*----------------------------   umfsplit.h     ---------------------------*/


#ifdef __WITH_UMFPACK__

#include  "iluinterface.h"
#include  "sparseblockmatrix.h"


#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#define HASHMAP   std::tr1::unordered_map
#define HASHSET   std::tr1::unordered_set
#else
#include  <ext/hash_map>
#include  <ext/hash_set>
#define HASHMAP  __gnu_cxx::hash_map
#define HASHSET  __gnu_cxx::hash_set
#endif
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
    class UmfSplit : virtual public IluInterface
  {
  private:
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
    UmfSplit(const MatrixInterface* A);
    UmfSplit();


    void SetMatrix(const MatrixInterface* A);
    
    ~UmfSplit();
    
    
    std::string GetName() const { return "Sparse-UMF"; }

    nvector<double>& GetRaw() { return __Ax; }
    vector<long>&    GetRawColumn() { return __Ac; }
    vector<long>&    GetRawPosition() { return __Ap; }
    
    
    int   n()          const { std::cerr << "UmfSplit::n()" << std::endl; abort(); }

    void zero()
    {
      __Ax.zero();
    }
    
    void ReInit(const SparseStructureInterface* SS)
    {
      assert(dynamic_cast<const SparseBlockMatrix<B> *>(__AS));
    }
    
    void copy_entries_fluid(const MatrixInterface*  A,
			    const vector<int>& fluid_l2g, 
			    const HASHSET<int>& interface_nodes);

    void copy_entries_solid(const MatrixInterface*  A,
			    const vector<int>& solid_l2g, 
			    const HASHSET<int>& interface_nodes);
    
    void copy_entries(const MatrixInterface*  A) {assert(0);};
    
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A) { assert(0);}
    void ConstructStructureFluid(const MatrixInterface& A, const vector<int>& fluid_l2g, const HASHSET<int>& interface_nodes);
    void ConstructStructureSolid(const MatrixInterface& A, const vector<int>& solid_l2g, const HASHSET<int>& interface_nodes);
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

#endif

/*----------------------------   umfsplit.h     ---------------------------*/
/* end of #ifndef __umfsplit_H */
#endif
/*----------------------------   umfsplit.h     ---------------------------*/
