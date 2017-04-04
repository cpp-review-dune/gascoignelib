/*----------------------------   aleblocks.h     ---------------------------*/
/*      $Id: aleblocks.h,v 1.5 2010/09/02 09:14:44 richter Exp $                 */
#ifndef __aleblocks_H
#define __aleblocks_H
/*----------------------------   aleblocks.h     ---------------------------*/

#include "fmatrixblock.h"
#include "numfixarray.h"


#define MATRIX_TYPE __MATRIX_TYPE__


    /**
     * Sparse Block for FSI-problems
     * only the necessary values are stored and used for calculation
     **/


namespace Gascoigne
{


  /**
   *  In Flow domain:
   *
   *   |                      |
   *   |  F[ns_ns]  FS[ns_u]  |
   *   |                      |
   *   |      0         S     |
   *   |                      |
   *
   * where A[ns_ns] in (DIM+1) x (DIM+1) matrix of the navier-stoeks system
   * and A[ns_u] in (DIM+1) x DIM derivative wrt deformation.
   * AD is laplace-matrix for extension of deformation. Just 1 value to be stored!
   *
   **/
  
  template<int DIM> 
    class FluidBlock
    {
      typedef nvector<double>::iterator        viterator;
      typedef nvector<double>::const_iterator  const_viterator;

      typedef typename numfixarray<(DIM+1)*(DIM+DIM+1),MATRIX_TYPE>::const_iterator  F_const_iterator;
      
  public:

      MATRIX_TYPE                                    __S;
      numfixarray<(DIM+1)*(DIM+DIM+1),MATRIX_TYPE>   __F;

      
      /////////// access
      int ncomp() const { return 2*DIM+1;}
      MATRIX_TYPE& diag(int i) 
	{
	  if (i<DIM+1) return __F[i*(DIM+DIM+1)+i];
	  return __S;
	}
      
      MATRIX_TYPE&       operator()(int i,int j)
	{
	  if (i<DIM+1) return __F[i*(DIM+DIM+1)+j];
	  if (i==j) return __S;
	  assert(0);
	}
      const MATRIX_TYPE& operator()(int i,int j) const
      {
	if (i<DIM+1) return __F[i*(DIM+DIM+1)+j];
	if (i==j) return __S;
	assert(0);
      }


      /////////// calculate
      void add(double s, const FluidBlock<DIM>& A)
      {
	__F.add(s,A.__F);
	__S += s * A.__S;
      }
      
      void add(double s, const TimePattern& TP){assert(0);}
      void zero()
      {
	__F.zero();
	__S=0;
      }
      void cadd(double s, viterator p, const_viterator q0) const
      {
	// p=p + s * (this) q

	//  p_v += s * F q_v
	viterator pp (p);
	F_const_iterator Fit = __F.begin();
	for (int i=0;i<DIM+1;++i,++pp)
	  {
	    const_viterator qq(q0);
	    for (int j=0;j<DIM+DIM+1;++j,++qq,++Fit)
	      *pp += s * *qq * (*Fit);
	  }
	
	// p_u += S q_u
	pp = p+DIM+1;
	const_viterator qq(q0+DIM+1);
	for (int i=0;i<DIM;++i,++pp,++qq)
	  *pp += s * __S * *qq;
      }
      
      void vector_get(nvector<MATRIX_TYPE>& v) const { assert(0); }
      void vector_set(nvector<MATRIX_TYPE>& v) { assert(0); }
      void vector_add(double d, nvector<MATRIX_TYPE>& v) { assert(0); }

      void   inverse ()
      {
	NodeMatrix<DIM+1,MATRIX_TYPE> X;
	for (int i=0;i<DIM+1;++i)
	  for (int j=0;j<DIM+1;++j)
	    X(i,j) = __F[i*(DIM+DIM+1)+j];
	X.inverse();
	for (int i=0;i<DIM+1;++i)
	  for (int j=0;j<DIM+1;++j)
	    __F[i*(DIM+DIM+1)+j] = X(i,j);
	assert(__S!=0);
	__S=1.0/__S;
	
	numfixarray< (DIM+1)*DIM , MATRIX_TYPE > Y;
	
	for (int i=0;i<DIM+1;++i)
	  for (int j=0;j<DIM;++j)
	    Y[i*DIM+j] = __F[i*(DIM+DIM+1)+j+DIM+1];
	for (int i=0;i<DIM+1;++i)
	  for (int j=0;j<DIM;++j)
	    {
	      __F[i*(DIM+DIM+1)+j+DIM+1] = 0;
	      for (int k=0;k<DIM+1;++k)
		__F[i*(DIM+DIM+1)+j+DIM+1] +=
		  - __S * __F[i*(DIM+DIM+1)+k] * Y[k*DIM+j];
	    }
      }
      inline void   vmult   (viterator p) const
      {
	numfixarray <DIM+DIM+1,double> vhelp;
	vhelp.zero();
	
	for (int i=0;i<DIM+DIM+1;++i)
	  vhelp[i] = *p++;
	p -= DIM+DIM+1;

	
	// p = F p_v + FS p_u
	F_const_iterator Fit = __F.begin();
	typename numfixarray<DIM+DIM+1,double>::const_iterator vhit = vhelp.begin();
	for (int i=0;i<DIM+1;++i,++p)
	  {
	    *p=0;

	    for (int j=0;j<DIM+DIM+1;++j,++vhit,++Fit)
	      *p += *Fit * *vhit;
	    vhit -= DIM+DIM+1;
	  }
	
	// += S p_u
	vhit += DIM+1;
	for (int i=0;i<DIM;++i,++p,++vhit)
	  *p += __S * *vhit;
      }

      FluidBlock<DIM>& operator+=(const FluidBlock<DIM>& v) { add(1.0,v); return *this; }
      FluidBlock<DIM>& operator*=(const FluidBlock<DIM>& v)
	{
	  numfixarray<(DIM+1)*(DIM+DIM+1),MATRIX_TYPE> save = __F;
	  __F.zero();
	  for (int i=0;i<DIM+1;++i)
	    for (int j=0;j<DIM+DIM+1;++j)
	      for (int k=0;k<DIM+1;++k)
		__F[i*(DIM+DIM+1) + j] += save[i*(DIM+DIM+1)+k] * v.__F[k*(DIM+DIM+1)+j];
	  for (int i=0;i<DIM+1;++i)
	    for (int j=0;j<DIM;++j)
	      __F[i*(DIM+DIM+1) + j+DIM+1] += __S * save[i*(DIM+DIM+1)+j+DIM+1];
	  __S *= v.__S;
	  return *this;
	}
      FluidBlock<DIM>& operator=(const FMatrixBlock<DIM+DIM+1>& X)
	{
	  for (int i=0;i<DIM+1;++i)
	    for (int j=0;j<DIM+DIM+1;++j)
	      __F[i*(DIM+DIM+1)+j] = X(i,j);
	  __S = X(DIM+1,DIM+1);
	  for (int j=0;j<DIM-1;++j)
	    assert(X(DIM+1,DIM+1)==X(DIM+1+1+j,DIM+1+1+j));
	  return *this;
	  
	}

      void subtract(viterator p0, const_viterator q0) const
      {
	cadd(-1.0,p0,q0);
      }
      
      void submult(const FluidBlock<DIM>& B, const FluidBlock<DIM>& C)
      {
	FluidBlock<DIM> X;
	X=B;
	X*=C;
	add(-1.0,X);
      }
      


      
      //////////// Boundary
      void   DirichletRow (const std::vector<int>& cv){assert(0);}
      void   DirichletCol (const std::vector<int>& cv){assert(0);}
      void   DirichletDiag(const std::vector<int>& cv){assert(0);}
      void   DirichletDiag(const std::vector<int>& cr,
			   const std::vector<int>& cc){assert(0);}

      void   getrow   (std::vector<double>& v, int i){assert(0);}
      void   getcolumn(std::vector<double>& v, int i){assert(0);}
      void   setrow   (std::vector<double>& v, int i){assert(0);}
      void   setcolumn(std::vector<double>& v, int i){assert(0);}


      //////////// entry
      void   entry     (const nmatrix<double>&){assert(0);}
      void   entry     (int i, int j, const EntryMatrix&, double s=1){assert(0);}
      void   dual_entry(int i, int j, const EntryMatrix&, double s=1.) { assert(0); }
      



      ////////// to be written when needed
      
      
      // just store, that matrix is transposed!
      void   transpose() { assert(0); }
      void   transpose(FluidBlock<DIM>& A) { assert(0); }
      
      friend std::ostream& operator<<(std::ostream &s,
				      const FluidBlock<DIM>& A)
	{ assert(0);}
      
      
      
    };


  //////////////////////////////////////////////////

  
  template<int DIM> 
    class SolidBlock
    {
      typedef nvector<double>::iterator        viterator;
      typedef nvector<double>::const_iterator  const_viterator;
      
  protected:
      FMatrixBlock<DIM+1>  __F;
      FMatrixBlock<DIM  >  __S;
      FMatrixBlock<DIM  >  __FS;

    public:

      /////////// access
      int ncomp() const { return 2*DIM+1;}
      MATRIX_TYPE& diag(int i) { 	assert(0); }
      MATRIX_TYPE&       operator()(int i,int j)       { assert(0);}
      const MATRIX_TYPE& operator()(int i,int j) const { assert(0);}


      /////////// calculate
      void add(double s, const FluidBlock<DIM>& A) { assert(0); }
      void add(double s, const TimePattern& TP){assert(0);}
      void zero(){assert(0);}
      void cadd(double s, viterator p, const_viterator q0) const { assert(0); }
      
      void vector_get(nvector<__MATRIX_TYPE__>& v) const { assert(0); }
      void vector_set(nvector<__MATRIX_TYPE__>& v) { assert(0); }
      void vector_add(double d, nvector<__MATRIX_TYPE__>& v) { assert(0); }

      void   inverse () { assert(0); }
      inline void   vmult   (viterator) const{assert(0);}

      FluidBlock<DIM>& operator+=(const FluidBlock<DIM>& v) { assert(0); }
      FluidBlock<DIM>& operator*=(const FluidBlock<DIM>& v) { assert(0); }

      void subtract(viterator p0, const_viterator q0) const { assert(0); }
      void   submult(const FluidBlock<DIM>& B, const FluidBlock<DIM>& C) { assert(0); }
      


      
      //////////// Boundary
      void   DirichletRow (const std::vector<int>& cv){assert(0);}
      void   DirichletCol (const std::vector<int>& cv){assert(0);}
      void   DirichletDiag(const std::vector<int>& cv){assert(0);}

      void   getrow   (std::vector<double>& v, int i){assert(0);}
      void   getcolumn(std::vector<double>& v, int i){assert(0);}
      void   setrow   (std::vector<double>& v, int i){assert(0);}
      void   setcolumn(std::vector<double>& v, int i){assert(0);}


      //////////// entry
      void   entry     (const nmatrix<double>&){assert(0);}
      void   entry     (int i, int j, const EntryMatrix&, double s=1){assert(0);}
      void   dual_entry(int i, int j, const EntryMatrix&, double s=1.) { assert(0); }
      



      ////////// to be written when needed
      
      
      // just store, that matrix is transposed!
      void   transpose() { assert(0); }
      void   transpose(FluidBlock<DIM>& A) { assert(0); }
      
      friend std::ostream& operator<<(std::ostream &s,
				      const FluidBlock<DIM>& A)
	{ assert(0);}
      
      
      
    };

    
  
  }


  
#include "sparseblockmatrix.h"
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

  template<class FD, class FS>
    void CopySubMatrixSolid(SparseBlockMatrix<FD >         * D,
			    const SparseBlockMatrix<FS > * S,
			    const std::vector<int>& nodes_l2g,
			    const HASHMAP<int,int>& nodes_g2l,
			    const HASHSET<int>& interface_nodes);
  template<class FD, class FS>
    void CopySubMatrixFluid(SparseBlockMatrix<FD >         * D,
			    const SparseBlockMatrix<FS > * S,
			    const std::vector<int>& nodes_l2g,
			    const HASHMAP<int,int>& nodes_g2l,
			    const HASHSET<int>& interface_nodes);
}



/*----------------------------   aleblocks.h     ---------------------------*/
/* end of #ifndef __aleblocks_H */
#endif
/*----------------------------   aleblocks.h     ---------------------------*/
