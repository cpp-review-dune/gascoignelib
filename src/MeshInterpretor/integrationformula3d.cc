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
