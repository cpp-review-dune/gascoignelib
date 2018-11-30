#include "umfsplit.h"
#include "fmatrixblock.h"
#include <fstream>
#ifdef __WITH_UMFPACK__

using namespace std;



#define UMFPACK_OK       0
#define UMFPACK_INFO    90
#define UMFPACK_CONTROL 20

#define UMFPACK_A	(0)	/* Ax=b		*/
#define UMFPACK_At	(1)	/* A'x=b	*/

namespace Gascoigne
{
  
  

extern "C" long umfpack_dl_symbolic
(
    long n,
    long m,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" long umfpack_dl_numeric
(
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" long umfpack_dl_solve
(
    long sys,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_dl_free_symbolic
(
    void **Symbolic
) ;
extern "C" void umfpack_dl_free_numeric
(
    void **Numeric
) ;

extern "C" long umfpack_dl_triplet_to_col
(
    long n,
    long nz,
    const long Ti [ ],
    const long Tj [ ],
    const double Tx [ ],
    long Bp [ ],
    long Bi [ ],
    double Bx [ ]
) ;

extern "C" void umfpack_dl_report_status
(
    const double Control [UMFPACK_CONTROL],
    long status
) ;
extern "C" void umfpack_dl_report_info
(
    const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_dl_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" long umfpack_dl_report_symbolic
(
    const char name [ ],
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" long umfpack_dl_report_numeric
(
    const char name [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_dl_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_dl_defaults
(
    const double Control [UMFPACK_CONTROL]
) ;

















// {
//   extern "C" int umfpack_di_symbolic
//   (
//    int n,
//    int m,
//    const int Ap [ ],
//    const int Ai [ ],
//    const double Ax [ ],
//    void **Symbolic,
//    const double Control [UMFPACK_CONTROL],
//    double Info [UMFPACK_INFO]
//    ) ;
//   extern "C" int umfpack_di_numeric
//   (
//    const int Ap [ ],
//    const int Ai [ ],
//    const double Ax [ ],
//    void *Symbolic,
//    void **Numeric,
//    const double Control [UMFPACK_CONTROL],
//    double Info [UMFPACK_INFO]
//    ) ;
//   extern "C" int umfpack_di_solve
//   (
//    int sys,
//    const int Ap [ ],
//    const int Ai [ ],
//    const double Ax [ ],
//    double X [ ],
//    const double B [ ],
//    void *Numeric,
//    const double Control [UMFPACK_CONTROL],
//    double Info [UMFPACK_INFO]
//    ) ;
//   extern "C" void umfpack_di_free_symbolic
//   (
//    void **Symbolic
//    ) ;
//   extern "C" void umfpack_di_free_numeric
//   (
//    void **Numeric
//    ) ;

//   extern "C" int umfpack_triplet_to_col
//   (
//    int n,
//    int nz,
//    const int Ti [ ],
//    const int Tj [ ],
//    const double Tx [ ],
//    int Bp [ ],
//    int Bi [ ],
//    double Bx [ ]
//    ) ;

//   extern "C" void umfpack_di_report_status
//   (
//    const double Control [UMFPACK_CONTROL],
//    int status
//    ) ;
//   extern "C" void umfpack_di_report_info
//   (
//    const double Control [UMFPACK_CONTROL],
//    const double Info [UMFPACK_INFO]
//    ) ;
//   extern "C" void umfpack_di_report_control
//   (
//    const double Control [UMFPACK_CONTROL]
//    ) ;
//   extern "C" int umfpack_di_report_symbolic
//   (
//    const char name [ ],
//    void *Symbolic,
//    const double Control [UMFPACK_CONTROL]
//    ) ;
//   extern "C" int umfpack_di_report_numeric
//   (
//    const char name [ ],
//    void *Numeric,
//    const double Control [UMFPACK_CONTROL]
//    ) ;
//   extern "C" void umfpack_di_report_control
//   (
//    const double Control [UMFPACK_CONTROL]
//    ) ;
//   extern "C" void umfpack_di_defaults
//   (
//    const double Control [UMFPACK_CONTROL]
//    ) ;


  //  ==================================================

  // This Constructor cannot work, if A is not supplied
  template<class B>
  UmfSplit<B>::UmfSplit()
    : __AS(0), Control(NULL), Info(NULL), 
      Symbolic(NULL), Numeric(NULL)
  {
    // __ncomp = 0;

    // Control = new double[UMFPACK_CONTROL];
    // umfpack_di_defaults(Control);
    // Control[0] = 2;
  }

  template<class B>
  UmfSplit<B>::UmfSplit(const MatrixInterface* A)
    : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
  {
    __AS = dynamic_cast<const SparseBlockMatrix<B>* >(A);
    assert(__AS);
    __ncomp = 0;

    Control = new double[UMFPACK_CONTROL];
    umfpack_dl_defaults(Control);
    Control[0] = 2;
  }
  template<class B>
  void UmfSplit<B>::SetMatrix(const MatrixInterface* A)
  {
    __AS = dynamic_cast<const SparseBlockMatrix<B>* >(A);
    assert(__AS);    

    __ncomp = 0;

    if (Control==0)
      {
	Control = new double[UMFPACK_CONTROL];
	umfpack_dl_defaults(Control);
	Control[0] = 2;
      }
  }
  
  
  template<class B>
  UmfSplit<B>::~UmfSplit()
  {
    umfpack_dl_free_symbolic (&Symbolic) ;
    umfpack_dl_free_numeric (&Numeric) ;
    
    if(Control) delete[] Control; Control=NULL;
    if(Info) delete[] Info; Info=NULL;
  }
  
  //  ==================================================


  template<class B>
  void UmfSplit<B>::ConstructStructureFluid(const MatrixInterface& A, 
					    const vector<int>& fluid_l2g, 
					    const HASHSET<int>& interface_nodes)
  {
    HASHSET<int> fluid_nodes;
    for (int i=0;i<fluid_l2g.size();++i)
      fluid_nodes.insert(fluid_l2g[i]);
    cout << "Fluid ILU: " << fluid_nodes.size() << " fluid nodes \t"
	 << interface_nodes.size() << " interface nodes" << endl;
    

    assert(__AS);
    // reinit size of UMFPACK-structure
    const ColumnStencil* ST = dynamic_cast<const ColumnStencil*> (__AS->GetStencil()); 
    // check if __AS and A are the same object...
    // then, we do not need __AS!!!
    assert(ST);
    __ncomp = __AS->mat(0)->ncomp();
    
    int n        = ST->n();
    int nentries = ST->nentries();

    // reserve size for matrix-entries
    __Ax.resize(__ncomp*__ncomp * nentries);
    __Ac.clear();
    __Ap.clear();
    
    __Ap.push_back(0);
    for (int r=0;r<n;++r)
      {
	for (int rc = 0 ; rc<__ncomp; ++rc)
	  {
	    for (int p = ST->start(r); p!=ST->stop(r); ++p)
	      {
		int c = ST->col(p);
		// insert column-values
		for (int cc = 0; cc<__ncomp;++cc)  __Ac.push_back(c*__ncomp+cc);
	      }
	    __Ap.push_back(__Ac.size());
	  }
      }

    long n_umf         = n * __ncomp;
    long n_entries_umf = nentries * __ncomp * __ncomp;

    assert(__Ac.size()==n_entries_umf);
    assert(__Ax.size()==n_entries_umf);
    assert(__Ap.size()==n_umf+1);
    assert(__Ap.size()>0);
    assert(__Ap[__Ap.size()-1] == n_entries_umf);

    umfpack_dl_free_symbolic (&Symbolic) ;
    
    const long* sb = &(*__Ap.begin());
    const long* cb = &(*__Ac.begin());
    
    int status = umfpack_dl_symbolic(n_umf, n_umf, sb,cb,
				     NULL, &Symbolic, Control, Info);
    
    if(status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control,Info);
      umfpack_dl_report_status(Control,status);
      ofstream file("MATRIX");
//       AP->Write(file);
      cerr << "umfpack_symbolic failed\n"; exit(1);
    }

  }


  template<class B>
  void UmfSplit<B>::ConstructStructureSolid(const MatrixInterface& A, const vector<int>& solid_l2g, const HASHSET<int>& interface_nodes) 
  {
    HASHSET<int> solid_nodes;
    for (int i=0;i<solid_l2g.size();++i)
      solid_nodes.insert(solid_l2g[i]);
    cout << "Solid ILU: " << solid_nodes.size() << " solid nodes \t"
	 << interface_nodes.size() << " interface nodes" << endl;
    
    assert(__AS);
    // reinit size of UMFPACK-structure
    const ColumnStencil* ST = dynamic_cast<const ColumnStencil*> (__AS->GetStencil()); 
    // check if __AS and A are the same object...
    // then, we do not need __AS!!!

    assert(ST);
    __ncomp = __AS->mat(0)->ncomp();
    
    int n        = ST->n();
    int nentries = ST->nentries();

    // reserve size for matrix-entries
    __Ax.resize(__ncomp*__ncomp * nentries);
    __Ac.clear();
    __Ap.clear();
    
    __Ap.push_back(0);
    for (int r=0;r<n;++r)
      {
	for (int rc = 0 ; rc<__ncomp; ++rc)
	  {
	    for (int p = ST->start(r); p!=ST->stop(r); ++p)
	      {
		int c = ST->col(p);
		// insert column-values
		for (int cc = 0; cc<__ncomp;++cc)  __Ac.push_back(c*__ncomp+cc);
	      }
	    __Ap.push_back(__Ac.size());
	  }
      }

    long n_umf         = n * __ncomp;
    long n_entries_umf = nentries * __ncomp * __ncomp;

    assert(__Ac.size()==n_entries_umf);
    assert(__Ax.size()==n_entries_umf);
    assert(__Ap.size()==n_umf+1);
    assert(__Ap.size()>0);
    assert(__Ap[__Ap.size()-1] == n_entries_umf);

    umfpack_dl_free_symbolic (&Symbolic) ;
    
    const long* sb = &(*__Ap.begin());
    const long* cb = &(*__Ac.begin());
    
    int status = umfpack_dl_symbolic(n_umf, n_umf, sb,cb,
				     NULL, &Symbolic, Control, Info);
    
    if(status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control,Info);
      umfpack_dl_report_status(Control,status);
      ofstream file("MATRIX");
//       AP->Write(file);
      cerr << "umfpack_symbolic failed\n"; exit(1);
    }

  }









  


  template<class B>
  void UmfSplit<B>::ConstructStructure(int ncomp, const SparseStructure& SS)
  {
    abort();
  }


  template<class B>
  void UmfSplit<B>::copy_entries_fluid(const MatrixInterface*  A,
				       const vector<int>& fluid_l2g, 
				       const HASHSET<int>& interface_nodes)
  {
    HASHSET<int> fluid_nodes;
    for (int i=0;i<fluid_l2g.size();++i)
      fluid_nodes.insert(fluid_l2g[i]);

    
	// check if __AS and A are the same object...
	// then, we do not need __AS!!!
	assert (static_cast<const void*> (__AS) == static_cast<const void*>  (A));
    
    // Copy Entries
    assert(__ncomp == __AS->mat(0)->ncomp());
    
    const ColumnStencil* ST = dynamic_cast<const ColumnStencil*> (__AS->GetStencil()); 
    assert(ST);
    int pp = 0;
    for (int r=0;r<ST->n();++r)
      {
	bool row_in_solid = (fluid_nodes.find(r)==fluid_nodes.end());
	bool row_on_inter = (interface_nodes.find(r)!=interface_nodes.end());
	
	for (int rc = 0 ; rc<__ncomp; ++rc)
	  {
	    for (int p = ST->start(r); p!=ST->stop(r); ++p)
	      {
		int col = ST->col(p);
		bool col_in_solid =  (fluid_nodes.find(col)==fluid_nodes.end());
		//		bool col_on_inter = (interface_nodes.find(col)!=interface_nodes.end());
		
		
		const B& b = *__AS->mat(p);
		for (int cc = 0; cc<__ncomp;++cc,++pp)
		  {
		    assert(pp<__Ax.size());
		    if (rc==0) // divergence equation
		      {
			if ((!row_in_solid)&&(!col_in_solid))
			  __Ax[pp] = b(rc,cc);
		      }
		    else // other equations. dirichlet on interface
		      {
			if ((!row_in_solid)&&(!col_in_solid))
			  if (!row_on_inter)
			    __Ax[pp] = b(rc,cc);
		      }
		    
		  }
	      }
	  }
      }
    
    // Set ones in diagonal of solid part 
    for (int r=0;r<__Ap.size()-1;++r)
      {
	int row  = r/__ncomp;
	//	int comp = r%__ncomp;
	// fluid -> do nothing
	if (fluid_nodes.find(row)!=fluid_nodes.end()) continue;
	

	// find diagonal
	for (int p=__Ap[r];p<__Ap[r+1];++p)
	  if (__Ac[p]==r)
	    {
	      __Ax[p]=1.0;
	    }
      }
    // Set some ones on interface
    for (int r=0;r<__Ap.size()-1;++r)
      {
	int row  = r/__ncomp;
	int comp = r%__ncomp;
	if (comp==0) continue;
	
	// fluid -> do nothing
	if (interface_nodes.find(row)==interface_nodes.end()) continue;
	
	
	// find diagonal
	for (int p=__Ap[r];p<__Ap[r+1];++p)
	  if (__Ac[p]==r)
	    {
	      __Ax[p]=1.0;
	    }
      }
  }


  template<class B>
  void UmfSplit<B>::copy_entries_solid(const MatrixInterface*  A,
				       const vector<int>& solid_l2g, 
				       const HASHSET<int>& interface_nodes)
    
  {
    HASHSET<int> solid_nodes;
    for (int i=0;i<solid_l2g.size();++i)
      solid_nodes.insert(solid_l2g[i]);
    
    // check if __AS and A are the same object...
    // then, we do not need __AS!!!
    assert (static_cast<const void*> (__AS) == static_cast<const void*>  (A));
    
    // Copy Entries
    assert(__ncomp == __AS->mat(0)->ncomp());
    
    const ColumnStencil* ST = dynamic_cast<const ColumnStencil*> (__AS->GetStencil()); 
    assert(ST);
    int pp = 0;
    for (int r=0;r<ST->n();++r)
      {
	bool row_in_fluid = (solid_nodes.find(r)==solid_nodes.end());

	for (int rc = 0 ; rc<__ncomp; ++rc)
	  {
	    for (int p = ST->start(r); p!=ST->stop(r); ++p)
	      {
		int col = ST->col(p);
		bool col_in_fluid =  (solid_nodes.find(col)==solid_nodes.end());
		
		const B& b = *__AS->mat(p);
		for (int cc = 0; cc<__ncomp;++cc,++pp)
		  {
		    assert(pp<__Ax.size());
		    if ((!row_in_fluid)&&(!col_in_fluid))
		    __Ax[pp] = b(rc,cc);
		  }
	      }
	  }
      }

    // Set ones in diagonal of fluid part
    for (int r=0;r<__Ap.size()-1;++r)
      {
	int row  = r/__ncomp;
	//	int comp = r%__ncomp;
	// solid -> do nothing
	if (solid_nodes.find(row)!=solid_nodes.end()) continue;
	

	// find diagonal
	for (int p=__Ap[r];p<__Ap[r+1];++p)
	  if (__Ac[p]==r)
	    {
	      __Ax[p]=1.0;
	    }
      }

  }

  // ==================================================

  template<class B>
  void UmfSplit<B>::compute_ilu ()
  {
    umfpack_dl_free_numeric (&Numeric) ;
    
    const long* sb = &(*__Ap.begin());
    const long* cb = &(*__Ac.begin());
    const double* mb = &(*__Ax.begin());
    int status = umfpack_dl_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info) ;
    if(status != UMFPACK_OK)
      {
	umfpack_dl_report_info(Control,Info);
	umfpack_dl_report_status(Control,status);
	cerr << "umfpack_numeric failed\n"; exit(1);
      }
  }
  
/*-----------------------------------------*/


  template<class B>
  void UmfSplit<B>::solve(GlobalVector& x) const 
  {
    assert(__ncomp == x.ncomp());
    assert(__Ap.size() == x.size()+1);

    
    const long* sb = &(*__Ap.begin());
    const long* cb = &(*__Ac.begin());
    const double* mb = &(*__Ax.begin());

    GlobalVector b = x;
    double* xb = &(*x.begin());
    const double* bb = &(*b.begin());
    int status = umfpack_dl_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info) ;
    
    if(status != UMFPACK_OK)
      {
	umfpack_dl_report_info(Control,Info);
	umfpack_dl_report_status(Control,status);
	cerr << "umfpack_dl_solve failed\n"; exit(1);
      }
  }
  

  template class UmfSplit<FMatrixBlock<1> >;
  template class UmfSplit<FMatrixBlock<2> >;
  template class UmfSplit<FMatrixBlock<3> >;
  template class UmfSplit<FMatrixBlock<4> >;
  template class UmfSplit<FMatrixBlock<5> >;
  template class UmfSplit<FMatrixBlock<7> >;

}


#endif
