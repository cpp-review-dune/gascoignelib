#include "galerkinintegratorq2.h"

namespace Gascoigne
{

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::GalerkinIntegratorQ2<DIM>() {}
  
/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::~GalerkinIntegratorQ2<DIM>() {}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegratorQ2<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new QuadGauss9;
      if (!MassFormulaPointer())     MassFormulaPointer() = new QuadGauss9; // ?????
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new QuadGauss16;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new LineGauss3;
    }
  else if (DIM==3)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new HexGauss27;
      if (!MassFormulaPointer())     MassFormulaPointer() = new HexGauss27; // ?????
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new HexGauss64;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new QuadGauss9;
    }
  assert(FormFormulaPointer());
  assert(ErrorFormulaPointer());
  assert(BoundaryFormulaPointer());
  assert(MassFormulaPointer());
}

/* ----------------------------------------- */

template GalerkinIntegratorQ2<2>;
template GalerkinIntegratorQ2<3>;

}
