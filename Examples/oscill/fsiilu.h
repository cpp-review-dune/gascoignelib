/*----------------------------   fsiilu.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsiilu_H
#define __fsiilu_H
/*----------------------------   fsiilu.h     ---------------------------*/





#include  "iluinterface.h"
#include  "splitbase.h"
#include  "sparse_umf.h"
#include  "sparseblockilu.h"
//#include  "sparseblockvanka.h"
#include  "fmatrixblock.h"
#include  "fsimatrix.h"
#include  "fsiilubase.h"
#include  "mginterpolatornested.h"
/*-------------------------------------------------------------*/



namespace Gascoigne
{

  template<int DIM>
    class FSIIlu : public IluInterface, public SplitBase, public FSIIluBase
  {  
  private:
    static constexpr int FCOMP=2*DIM+1;
    
    const FSIMatrix<DIM>* __FSIMATRIX;
    
    //SparseUmf<FMatrixBlock<DIM> > __SM;
    SparseBlockIlu<DiagBlock<DIM> > __SM;
    SparseUmf<FMatrixBlock<DIM> > __SES;
    //    SparseBlockIlu<FMatrixBlock<DIM> > __SES;
    //SparseBlockVanka<DIM,DIM> __SES;
    
     
    //SparseUmf<FMatrixBlock<DIM+1> >  __F_NS;
    SparseBlockIlu<FMatrixBlock<DIM+1> > __F_NS;
      
    SparseBlockIlu<DiagBlock<DIM > >  __F_EXT;
    //    SparseDiagUmf<DIM> __F_EXT;
      

    mutable GlobalVector _xs1,_ys1,_hs1, _xs2,_ys2,_hs2;
    mutable GlobalVector _xf_ns,_yf_ns, _hf_ns, _xf_ext, _yf_ext;
    mutable GlobalVector _hf_ext;






    
  public:
    // For multigrid
    int thislevel;
    int coarselevel;
    vector<const FSIIlu<DIM> *> __FSI_H;
    MgInterpolatorNested     __interpolator;
    nvector<int> __boundarynodes;
    

    void init_mg(int l,int cl, vector<const FSIIlu<DIM> *> ILUH, const MgInterpolatorNested& mgint, const nvector<int>& bn)
    {
      //      assert(bn.size()==0);
      thislevel   = l;
      coarselevel = cl;
      __FSI_H.resize(ILUH.size());
      for (int i=0;i<ILUH.size();++i) __FSI_H[i]=ILUH[i];
      __interpolator = mgint;
      __boundarynodes = bn;
      /* cout << l << " " << cl << "\t"; */
      /* cout << __interpolator.GetZweier().size() << "\t"; */
      /* cout << __interpolator.GetVierer().size() << "\t"; */
      /* cout << __interpolator.GetAchter().size() << "\t"; */
      /* cout << __interpolator.GetC2F().size() << endl; */
      
    }
    


    void mg_smooth_ses(GlobalVector& x,GlobalVector& y, GlobalVector& h, int iter) const;
    void mg_ses(GlobalVector& x,GlobalVector& y, GlobalVector& h) const;
    void mg_restrict_ses(GlobalVector& Y,const GlobalVector& h) const;
    void mg_prolongate_add_ses(GlobalVector& Y,const GlobalVector& h) const;
    
      
    FSIIlu() {}
    FSIIlu(const MatrixInterface* M);
      
    virtual ~FSIIlu() {};

    //const SparseUmf<FMatrixBlock<DIM> >& GetSM() const { return __SM; }
    const SparseBlockIlu<DiagBlock<DIM> >& GetSM() const { return __SM; }
    const SparseUmf<FMatrixBlock<DIM> >& GetSES() const { return __SES; }
    //    const SparseBlockVanka<DIM,DIM>& GetSES() const { return __SES; }
    //const SparseBlockIlu<FMatrixBlock<DIM> >& GetSES() const { return __SES; }

      
    const SparseBlockIlu<FMatrixBlock<DIM+1> >& GetF_NS()  const { return __F_NS; }
    //    const SparseUmf<FMatrixBlock<DIM+1> >& GetF_NS()  const { return __F_NS; }
    
    const SparseBlockIlu<DiagBlock<DIM> >& GetF_EXT() const { return __F_EXT; }
    //    const SparseDiagUmf<DIM >& GetF_EXT() const { return __F_EXT; }


    void modify_ilu_F_NS()
    {
      __F_NS.modify(0,0.01);
      for (int c=1;c<=DIM;++c)
	__F_NS.modify(c,0.1);
    }
    void modify_ilu_S()
    {
      /* for (int c=0;c<DIM;++c) */
      /* 	__SM.modify(c,0.1); */
      /* for (int c=0;c<DIM;++c)  */
      /* 	__SES.modify(c,0.); */
    }
      
    int   n() const{assert(0);}
    std::string GetName() const{assert(0);}
    void ReInit(const SparseStructureInterface* A);
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A){std::cerr << "!CS!" << std::endl; abort();}
    void ConstructStructure(const IntVector& permlfuid, const IntVector& permsolid, const MatrixInterface& A);

      
    /* void Fluid_g2l_set(GlobalVector& xf, const GlobalVector& x) const; */
    void Solid1_g2l_set(GlobalVector& xs, const GlobalVector& x) const;
    void Solid2_g2l_set(GlobalVector& xs, const GlobalVector& x) const;
    /* void Fluid_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const; */
    void Solid1_l2g_add(GlobalVector& x, const GlobalVector& xs,double s=1.0) const;
    void Solid2_l2g_add(GlobalVector& x, const GlobalVector& xs,double s=1.0) const;
      
    void Fluid_NS_g2l_set(GlobalVector& xf, const GlobalVector& x) const;
    void Fluid_EXT_g2l_set(GlobalVector& xf, const GlobalVector& x) const;
    void Fluid_NS_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const;
    void Fluid_EXT_l2g_add(GlobalVector& x, const GlobalVector& xf,double s=1.0) const;

      
    void zero()
    {
      //	__F.zero();
      __SM.zero();
      __SES.zero();

      __F_NS.zero();
      __F_EXT.zero();
    }
    void compute_ilu ();
    void copy_entries(const MatrixInterface* A);

      
    void Precondition(const IluInterface& M,GlobalVector &x) const;


    // vmult that does the correct interface handling even tho A has no dirichlet values embedded!
    void  vmult_NS_dirichlet(const MatrixInterface &A, 
			     GlobalVector &y, GlobalVector &x, double s) const ;
    void  vmult_EXT_dirichlet(const MatrixInterface &A, 
			     GlobalVector &y, GlobalVector &x, double s) const ;
      
    void  Jacobi(const MatrixInterface &A, 
		 GlobalVector &x) const;
      

    /* int  BiCGSTABL(int L, */
    /* 		   const MatrixInterface &A,  */
    /* 		   GlobalVector &x, const GlobalVector &b, */
    /* 		   const IluInterface &M,  */
    /* 		   int &max_iter, double &tol) const ; */

    int  BiCGSTAB(const MatrixInterface &A, 
		  GlobalVector &x, const GlobalVector &b,
		  const IluInterface &M, 
		  int &max_iter, double &tol) const ;
    int  BiCGSTAB(const MatrixInterface &A1, const MatrixInterface &A2, 
		  GlobalVector &x, const GlobalVector &b,
		  const IluInterface &M, 
		  int &max_iter, double &tol) const ;
    int  CG(const MatrixInterface &A, 
	    GlobalVector &x, const GlobalVector &b,
	    const IluInterface &M, 
	    int &max_iter, double &tol) const ;
      
    int  CGSolid(const MatrixInterface &A, 
	    GlobalVector &x, const GlobalVector &b,
	    const IluInterface &M, 
	    int &max_iter, double &tol) const ;
      

      
      
    void solve(GlobalVector& x) const;
    void solve_fluid(GlobalVector& x) const;
    void solve_solid(GlobalVector& x) const;
      
    void solve_transpose(GlobalVector& x) const {
      std::cerr << "\"IluInterface::solve_transpose\" not written!" << std::endl;
      abort();
    }
    std::ostream& Write(std::ostream &s) const {
      std::cerr << "\"IluInterface::Write\" not written!" << std::endl;
      abort();
    }
  };
}






/*----------------------------   fsiilu.h     ---------------------------*/
/* end of #ifndef __fsiilu_H */
#endif
/*----------------------------   fsiilu.h     ---------------------------*/
