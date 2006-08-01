#ifndef  __GalerkinIntegratorQ4_h
#define  __GalerkinIntegratorQ4_h

#include  "galerkinintegrator.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinIntegratorQ4

////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinIntegratorQ4 : virtual public GalerkinIntegrator<DIM>
{
public:

//
////  Con(De)structor 
//
  GalerkinIntegratorQ4<DIM>();

    ~GalerkinIntegratorQ4<DIM>();
  std::string GetName() const {return "GalerkinIntegratorQ4";}

  void BasicInit();
};
}

#endif
