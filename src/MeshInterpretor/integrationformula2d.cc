#include  "integrationformula.h"
#include  "integrationformula.xx"
#include  "tensorformula2d.xx"

/*------------------------------------------------------------*/

template IntegrationFormula<1>;
template IntegrationFormula<2>;

/*------------------------------------------------------------*/

template TensorFormula2d<1,LineMidPoint>;
template TensorFormula2d<2,LineTrapez  >;
template TensorFormula2d<3,LineSimpson >;
template TensorFormula2d<1,LineGauss1  >;
template TensorFormula2d<2,LineGauss2  >;
template TensorFormula2d<3,LineGauss3  >;
template TensorFormula2d<4,LineGauss4  >;
template TensorFormula2d<5,LineGauss5  >;
template TensorFormula2d<6,LineGauss6  >;
template TensorFormula2d<7,LineGauss7  >;
template TensorFormula2d<8,LineGauss8  >;
template TensorFormula2d<9,LineGauss9  >;
template TensorFormula2d<10,LineGauss10>;

/*------------------------------------------------------------*/
