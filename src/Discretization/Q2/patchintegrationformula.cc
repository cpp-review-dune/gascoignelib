#include  "patchintegrationformula.h"
#include  "patchintegrationformula.xx"

namespace Gascoigne
{

/*------------------------------------------------------------*/

template class PatchFormula1d<3,LineGauss3>;
template class PatchFormula1d<4,LineGauss4>;

/*------------------------------------------------------------*/

template class PatchFormula2d<4,TensorFormula2d<2,LineTrapez> >;
template class PatchFormula2d<16,TensorFormula2d<4,LineGauss4> >;
template class PatchFormula2d<9,TensorFormula2d<3,LineGauss3> >;
template class PatchFormula2d<4,TensorFormula2d<2,LineGauss2> >;

/*------------------------------------------------------------*/

template class PatchFormula3d<8, TensorFormula3d<2,LineGauss2> >;
template class PatchFormula3d<27,TensorFormula3d<3,LineGauss3> >;
template class PatchFormula3d<64,TensorFormula3d<4,LineGauss4> >;
template class PatchFormula3d<125,TensorFormula3d<5,LineGauss5> >;

/*------------------------------------------------------------*/

}

