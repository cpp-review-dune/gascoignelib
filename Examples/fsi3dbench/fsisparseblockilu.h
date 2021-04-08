/*----------------------------   fsisparseblockilu.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __fsisparseblockilu_H
#define __fsisparseblockilu_H
/*----------------------------   fsisparseblockilu.h
 * ---------------------------*/

#include "gascoignehash.h"
#include "iluinterface.h"
#include "sparseblockilu.h"
/*-------------------------------------------------------------*/

namespace Gascoigne {
template <class B> class FSISparseBlockIlu : public virtual SparseBlockIlu<B> {
protected:
public:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  void copy_entries(const HASHMAP<int, int> &G2L, const HASHSET<int> &INT,
                    const MatrixInterface *A);
#pragma GCC diagnostic pop
};
} // namespace Gascoigne

/*----------------------------   fsisparseblockilu.h
 * ---------------------------*/
/* end of #ifndef __fsisparseblockilu_H */
#endif
/*----------------------------   fsisparseblockilu.h
 * ---------------------------*/
