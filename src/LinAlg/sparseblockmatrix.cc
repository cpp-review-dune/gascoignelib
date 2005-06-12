#include "sparseblockmatrix.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template class SparseBlockMatrix<FMatrixBlock<1> >;
template class SparseBlockMatrix<FMatrixBlock<2> >;
template class SparseBlockMatrix<FMatrixBlock<3> >;
template class SparseBlockMatrix<FMatrixBlock<4> >;


template class SparseBlockMatrix<CFDBlock3d>;
}
