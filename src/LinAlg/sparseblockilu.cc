#include  "sparseblockilu.xx"
#include  "fmatrixblock.h"
#include  "cfdblock3d.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template SparseBlockIlu<FMatrixBlock<1> >;
template SparseBlockIlu<FMatrixBlock<2> >;
template SparseBlockIlu<FMatrixBlock<3> >;
template SparseBlockIlu<FMatrixBlock<4> >;
template SparseBlockIlu<FMatrixBlock<5> >;
template SparseBlockIlu<FMatrixBlock<6> >;
template SparseBlockIlu<FMatrixBlock<7> >;
template SparseBlockIlu<FMatrixBlock<8> >;
template SparseBlockIlu<FMatrixBlock<9> >;
template SparseBlockIlu<FMatrixBlock<10> >;
template SparseBlockIlu<FMatrixBlock<11> >;
template SparseBlockIlu<FMatrixBlock<12> >;


template SparseBlockIlu<CFDBlock3d>;
}
