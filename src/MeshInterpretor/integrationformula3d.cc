#include  "integrationformula.h"
#include  "integrationformula.xx"
#include  "tensorformula3d.xx"

/*------------------------------------------------------------*/

template class IntegrationFormula<3>;

template class TensorFormula3d<2,LineTrapez>;
template class TensorFormula3d<2,LineGauss2>;
template class TensorFormula3d<3,LineGauss3>;
template class TensorFormula3d<4,LineGauss4>;

/*------------------------------------------------------------*/
