#include  "sparseblockilu.xx"
#include  "fmatrixblock.h"
#include  "cfdblock3d.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template class SparseBlockIlu<FMatrixBlock<1> >;
template class SparseBlockIlu<FMatrixBlock<2> >;
template class SparseBlockIlu<FMatrixBlock<3> >;
template class SparseBlockIlu<FMatrixBlock<4> >;
template class SparseBlockIlu<FMatrixBlock<5> >;
template class SparseBlockIlu<FMatrixBlock<6> >;


template class SparseBlockIlu<CFDBlock3d>;
}
