#include  "integrationformula.h"
#include  "integrationformula.xx"
#include  "tensorformula2d.xx"

/*------------------------------------------------------------*/

template class IntegrationFormula<1>;
template class IntegrationFormula<2>;

/*------------------------------------------------------------*/

template class TensorFormula2d<1,LineMidPoint>;
template class TensorFormula2d<2,LineTrapez  >;
template class TensorFormula2d<3,LineSimpson >;
template class TensorFormula2d<1,LineGauss1  >;
template class TensorFormula2d<2,LineGauss2  >;
template class TensorFormula2d<3,LineGauss3  >;
template class TensorFormula2d<4,LineGauss4  >;
template class TensorFormula2d<5,LineGauss5  >;
template class TensorFormula2d<6,LineGauss6  >;
template class TensorFormula2d<7,LineGauss7  >;
template class TensorFormula2d<8,LineGauss8  >;
template class TensorFormula2d<9,LineGauss9  >;
template class TensorFormula2d<10,LineGauss10>;

/*------------------------------------------------------------*/
