#include "galerkinintegratorq4.h"

namespace Gascoigne
{

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ4<DIM>::GalerkinIntegratorQ4() {}

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ4<DIM>::~GalerkinIntegratorQ4<DIM>() {}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegratorQ4<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new QuadGauss16;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new QuadGauss25; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new QuadGauss25;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new LineGauss4;
    }
  else if (DIM==3)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new HexGauss64;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new HexGauss125; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new HexGauss125;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new QuadGauss16;
    }
  assert(GalerkinIntegrator<DIM>::FormFormulaPointer());
  assert(GalerkinIntegrator<DIM>::ErrorFormulaPointer());
  assert(GalerkinIntegrator<DIM>::BoundaryFormulaPointer());
  assert(GalerkinIntegrator<DIM>::MassFormulaPointer());
}

/* ----------------------------------------- */

template class GalerkinIntegratorQ4<2>;
template class GalerkinIntegratorQ4<3>;

}
