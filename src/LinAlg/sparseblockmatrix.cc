#include "sparseblockmatrix.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template SparseBlockMatrix<FMatrixBlock<1> >;
template SparseBlockMatrix<FMatrixBlock<2> >;
template SparseBlockMatrix<FMatrixBlock<3> >;
template SparseBlockMatrix<FMatrixBlock<4> >;
template SparseBlockMatrix<FMatrixBlock<5> >;
template SparseBlockMatrix<FMatrixBlock<6> >;

template SparseBlockMatrix<CFDBlock3d>;
}
