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
  GalerkinIntegratorQ2<DIM>();
  
    ~GalerkinIntegratorQ2<DIM>();
  std::string GetName() const {return "GalerkinIntegratorQ2";}

  void BasicInit();
};
}

#endif
