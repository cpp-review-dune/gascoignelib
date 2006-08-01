#include  "integrationformula.h"
#include  "integrationformulasummed.h"
#include  "tensorformula2d.xx"
#include  "tensorformula3d.xx"

/*------------------------------------------------------------*/

namespace Gascoigne
{
template class IntegrationFormulaBase<1>;
template class IntegrationFormulaBase<2>;
template class IntegrationFormulaBase<3>;

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

template class TensorFormula3d<1,LineGauss1>;
template class TensorFormula3d<2,LineTrapez>;
template class TensorFormula3d<2,LineGauss2>;
template class TensorFormula3d<3,LineGauss3>;
template class TensorFormula3d<4,LineGauss4>;
template class TensorFormula3d<5,LineGauss5>;

/*------------------------------------------------------------*/
template class IntegrationFormulaSummed2d<QuadGauss1>;
template class IntegrationFormulaSummed2d<QuadGauss4>;
template class IntegrationFormulaSummed3d<QuadGauss1>;
}

