#include  "threadsparseblockilu.xx"
#include  "fmatrixblock.h"
#include  "cfdblock3d.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template class ThreadSparseBlockIlu<FMatrixBlock<1> >;
template class ThreadSparseBlockIlu<FMatrixBlock<2> >;
template class ThreadSparseBlockIlu<FMatrixBlock<3> >;
template class ThreadSparseBlockIlu<FMatrixBlock<4> >;
template class ThreadSparseBlockIlu<FMatrixBlock<5> >;
template class ThreadSparseBlockIlu<FMatrixBlock<6> >;
template class ThreadSparseBlockIlu<FMatrixBlock<7> >;
template class ThreadSparseBlockIlu<FMatrixBlock<8> >;
}
