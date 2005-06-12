#include "galerkinintegratorq2.h"

namespace Gascoigne
{

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::GalerkinIntegratorQ2() {}
  
/* ----------------------------------------- */

template<int DIM>
GalerkinIntegratorQ2<DIM>::~GalerkinIntegratorQ2<DIM>() {}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegratorQ2<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new QuadGauss9;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new QuadGauss9; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new QuadGauss16;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new LineGauss3;
    }
  else if (DIM==3)
    {
      if (!GalerkinIntegrator<DIM>::FormFormulaPointer())     GalerkinIntegrator<DIM>::FormFormulaPointer() = new HexGauss27;
      if (!GalerkinIntegrator<DIM>::MassFormulaPointer())     GalerkinIntegrator<DIM>::MassFormulaPointer() = new HexGauss27; // ?????
      if (!GalerkinIntegrator<DIM>::ErrorFormulaPointer())    GalerkinIntegrator<DIM>::ErrorFormulaPointer() = new HexGauss64;
      if (!GalerkinIntegrator<DIM>::BoundaryFormulaPointer()) GalerkinIntegrator<DIM>::BoundaryFormulaPointer() = new QuadGauss9;
    }
  assert(GalerkinIntegrator<DIM>::FormFormulaPointer());
  assert(GalerkinIntegrator<DIM>::ErrorFormulaPointer());
  assert(GalerkinIntegrator<DIM>::BoundaryFormulaPointer());
  assert(GalerkinIntegrator<DIM>::MassFormulaPointer());
}

/* ----------------------------------------- */

template class GalerkinIntegratorQ2<2>;
template class GalerkinIntegratorQ2<3>;

}
