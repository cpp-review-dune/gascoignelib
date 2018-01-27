#include "fsisparse_umf.h"
#include "fmatrixblock.h"

namespace Gascoigne
{
  
  template<class B>
  void FSISparseUmf<B>::copy_entries(const HASHMAP<int,int>& G2L,
				     const HASHSET<int>& INT,
				     const MatrixInterface*  A)
  {
    // check if __AS and A are the same object...
    // then, we do not need __AS!!!
    assert (static_cast<const void*> (this->__AS) == static_cast<const void*>  (A));
    
    // Copy Entries
    assert(this->__ncomp == this->__AS->mat(0)->ncomp());
    
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (this->__AS->GetStencil()); 
    assert(ST);
    int pp = 0;
    for (int r=0;r<ST->n();++r)
      {
	// only rows in g2l
	auto it_row = G2L.find(r);
	bool row_found = (it_row!=G2L.end());
	bool on_int    = (INT.find(r)!=INT.end());
	
	for (int rc = 0 ; rc<this->__ncomp; ++rc)
	  {
	    for (int p = ST->start(r); p!=ST->stop(r); ++p)
	      {
		int col = ST->col(p);
		auto it_col = G2L.find(col);
		bool col_found = (it_col!=G2L.end());

		const B& b = *this->__AS->mat(p);
		for (int cc = 0; cc<this->__ncomp;++cc,++pp)
		  {
		    // all fluid-fluid but not interface-fluid for velocity
		    if ((col_found)&&(row_found) && !( (rc>0)&&(on_int) ))
		      {
			assert(pp<this->__Ax.size());
			this->__Ax[pp] = b(rc,cc);
		      }
		    else this->__Ax[pp]=0.0;
		    
		    if (!(row_found) || ((rc>0) && (on_int) ) )
		      {
			if (r==col)
			  if (rc==cc)
			    this->__Ax[pp]=1.0;
		      }
		  }
	      }
	  } 
      }
    
  } 
  
  
  template class FSISparseUmf<FMatrixBlock<3> >;
  template class FSISparseUmf<FMatrixBlock<4> >;
  
    
}
