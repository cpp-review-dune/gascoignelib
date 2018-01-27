/*----------------------------   fsisparse_umf.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsisparse_umf_H
#define __fsisparse_umf_H
/*----------------------------   fsisparse_umf.h     ---------------------------*/

#include "sparse_umf.h"
#include "gascoignehash.h"
namespace Gascoigne
{
  template<class B>
  class FSISparseUmf : virtual public SparseUmf<B>
  {
  public:
    void copy_entries(const HASHMAP<int,int>& G2L,
		      const HASHSET<int>& INT,
		      const MatrixInterface*  A);
    
  };

}



/*----------------------------   fsisparse_umf.h     ---------------------------*/
/* end of #ifndef __fsisparse_umf_H */
#endif
/*----------------------------   fsisparse_umf.h     ---------------------------*/
