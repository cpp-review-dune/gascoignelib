#include  "integrationformula.h"
#include  "integrationformula.xx"
#include  "tensorformula3d.xx"

/*------------------------------------------------------------*/

template IntegrationFormula<3>;

template TensorFormula3d<2,LineTrapez>;
template TensorFormula3d<2,LineGauss2>;
template TensorFormula3d<3,LineGauss3>;
template TensorFormula3d<4,LineGauss4>;

/*------------------------------------------------------------*/

template PatchFormula2d<9,TensorFormula2d<3,LineGauss3> >;
template PatchFormula2d<4,TensorFormula2d<2,LineGauss2> >;
template PatchFormula3d<8, TensorFormula3d<2,LineGauss2> >;
template PatchFormula3d<27,TensorFormula3d<3,LineGauss3> >;

/*------------------------------------------------------------*/
