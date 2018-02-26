/*----------------------------   fsimatrix.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsimatrix_H
#define __fsimatrix_H
/*----------------------------   fsimatrix.h     ---------------------------*/


#include  "matrixinterface.h"
#include  "fmatrixblock.h"
#include  "fmatrixblock_nm.h"
#include  "diagblock.h"
#include  "sparseblockmatrix.h"
#include  "splitbase.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

  template<int DIM>
    class FSIMatrix : public MatrixInterface, public SplitBase
  {
    static constexpr int FCOMP=2*DIM+1;
    
    //    SparseBlockMatrix<FMatrixBlock<FCOMP> > __F;
    SparseBlockMatrix<DiagBlock<DIM> > __S11,__S21;
    SparseBlockMatrix<FMatrixBlock<DIM> > __S12,__S22;

    SparseBlockMatrix<FMatrixBlock<DIM+1> >        __F_NS;
    SparseBlockMatrix<DiagBlock<DIM> >             __F_EXT;
    SparseBlockMatrix<FMatrixBlockNM<DIM+1,DIM> >  __F_ALE;

    mutable GlobalVector _xs1,_xs2,_ys1,_ys2;

    mutable GlobalVector _xf_ns,_yf_ns, _xf_ext, _yf_ext;

  public:

    
    FSIMatrix() { }
    virtual ~FSIMatrix() { }

    virtual std::string GetName() const { return "FSI Matrix"; }

    //    const SparseBlockMatrix<FMatrixBlock<FCOMP> >& GetF() const { return __F; }
    const SparseBlockMatrix<DiagBlock<DIM> >& GetS11() const { return __S11; }
    const SparseBlockMatrix<FMatrixBlock<DIM> >& GetS12() const { return __S12; }
    const SparseBlockMatrix<DiagBlock<DIM> >& GetS21() const { return __S21; }
    const SparseBlockMatrix<FMatrixBlock<DIM> >& GetS22() const { return __S22; }

    const SparseBlockMatrix<FMatrixBlock<DIM+1> >      & GetF_NS()  const { return __F_NS; }
    const SparseBlockMatrix<DiagBlock<DIM  > >         & GetF_EXT() const { return __F_EXT; }
    const SparseBlockMatrix<FMatrixBlockNM<DIM+1,DIM> >& GetF_ALE() const { return __F_ALE; }

    const MatrixInterface* GetFluid_NS()  const { return &__F_NS; }
    const MatrixInterface* GetFluid_EXT() const { return &__F_EXT; }
    const MatrixInterface* GetFluid_ALE() const { return &__F_ALE; }
    const MatrixInterface* GetSolid11()     const { return &__S11; }
    const MatrixInterface* GetSolid12()     const { return &__S12; }
    const MatrixInterface* GetSolid21()     const { return &__S21; }
    const MatrixInterface* GetSolid22()     const { return &__S22; }
    
    
    void Fluid_g2l_set(GlobalVector& xf, const GlobalVector& x) const;
    void Solid1_g2l_set(GlobalVector& xs, const GlobalVector& x) const;
    void Solid2_g2l_set(GlobalVector& xs, const GlobalVector& x) const;
    void Fluid_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const;
    void Solid1_l2g_add(GlobalVector& x, const GlobalVector& xs,double s=1.0) const;
    void Solid2_l2g_add(GlobalVector& x, const GlobalVector& xs,double s=1.0) const;

    void Fluid_NS_g2l_set(GlobalVector& xf, const GlobalVector& x) const;
    void Fluid_EXT_g2l_set(GlobalVector& xf, const GlobalVector& x) const;
    void Fluid_NS_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const;
    void Fluid_EXT_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const;



    
    virtual const StencilInterface* GetStencil() const {assert(0);}
    virtual void ReInit(const SparseStructureInterface* S);
    

    void copy_entries(const MatrixInterface& S)
    {
      std::cerr << "\"MatrixInterface::copy_entries\" not written!" << std::endl;
      abort();
    }
    void zero()
    {
      //      __F.zero();
      __S11.zero();
      __S12.zero();
      __S21.zero();
      __S22.zero();
      __F_NS.zero();
      __F_EXT.zero();
      __F_ALE.zero();
    }
    void transpose()   { abort(); }
    std::ostream& Write(std::ostream& os) const { abort(); }

    //
    /// for matrix assembling
    //
    typedef IntVector::const_iterator niiterator;
    
    void entry(nvector<int>::const_iterator start1,
	       nvector<int>::const_iterator stop1,
	       nvector<int>::const_iterator start2,
	       nvector<int>::const_iterator stop2,
	       const EntryMatrix& M, double s=1.)
    { 
      std::cerr << "\"MatrixInterface::entry\" not written!" << std::endl;
      abort();
    }
    void entry(niiterator start, niiterator stop,
	       const EntryMatrix& M, double s=1.);
    
    void entry_diag(int i, const nmatrix<double>& M)
    { abort();}
      
    void vmult(GlobalVector& y, const GlobalVector& x, double s=1.) const;
    

    ////////// BOUNDARY
    void dirichlet (const nvector<int>& bv, const std::vector<int>& cv);
    void dirichlet (int i, const std::vector<int>& cv);
    

  };



  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////


  
  
}




/*----------------------------   fsimatrix.h     ---------------------------*/
/* end of #ifndef __fsimatrix_H */
#endif
/*----------------------------   fsimatrix.h     ---------------------------*/
