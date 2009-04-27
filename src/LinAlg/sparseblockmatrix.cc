#include "sparseblockmatrix.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template class SparseBlockMatrix<FMatrixBlock<1> >;
template class SparseBlockMatrix<FMatrixBlock<2> >;
template class SparseBlockMatrix<FMatrixBlock<3> >;
template class SparseBlockMatrix<FMatrixBlock<4> >;
template class SparseBlockMatrix<FMatrixBlock<5> >;
template class SparseBlockMatrix<FMatrixBlock<6> >;
template class SparseBlockMatrix<FMatrixBlock<7> >;
template class SparseBlockMatrix<FMatrixBlock<8> >;


template class SparseBlockMatrix<CFDBlock3d>;
}
