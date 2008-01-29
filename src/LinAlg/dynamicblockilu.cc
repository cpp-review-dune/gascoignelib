#include "dynamicblockilu.xx"
#include "fmatrixblock.h"
#include "cfdblock3d.h"

namespace Gascoigne
{
template class DynamicBlockIlu<FMatrixBlock<1> >;
template class DynamicBlockIlu<FMatrixBlock<2> >;
template class DynamicBlockIlu<FMatrixBlock<3> >;
template class DynamicBlockIlu<FMatrixBlock<4> >;


template class DynamicBlockIlu<CFDBlock3d>;
}
