#include "dynamicblockmatrix.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template class DynamicBlockMatrix<FMatrixBlock<1> >;
template class DynamicBlockMatrix<FMatrixBlock<2> >;
template class DynamicBlockMatrix<FMatrixBlock<3> >;
template class DynamicBlockMatrix<FMatrixBlock<4> >;


template class DynamicBlockMatrix<CFDBlock3d>;
}
