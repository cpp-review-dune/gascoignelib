#ifndef  __GalerkinIntegratorQ2_h
#define  __GalerkinIntegratorQ2_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinIntegratorQ2

////
////
/////////////////////////////////////////////

#include  "galerkinintegrator.h"

namespace Gascoigne
{

template<int DIM>
class GalerkinIntegratorQ2 : virtual public GalerkinIntegrator<DIM>
{
public:

//
////  Con(De)structor 
//
  GalerkinIntegratorQ2<DIM>() {
  if (DIM==2)
    {
      FormFormulaPointer() = new QuadGauss9;
      ErrorFormulaPointer() = new QuadGauss16;
      BoundaryFormulaPointer() = new LineGauss3;
    }
  else if (DIM==3)
    {
      FormFormulaPointer() = new HexGauss27;
      ErrorFormulaPointer() = new HexGauss64;
      BoundaryFormulaPointer() = new QuadGauss9;
    }
  assert(FormFormulaPointer());
  assert(ErrorFormulaPointer());
  assert(BoundaryFormulaPointer());
  }
  
  ~GalerkinIntegratorQ2<DIM>() {}
};
}

#endif
