#ifndef __IntegratorLpsWithSecond_h
#define __IntegratorLpsWithSecond_h

#include "integratorwithsecond.h"
#include "galerkinlpsintegratorq2.h"

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
class IntegratorLpsWithSecond : public IntegratorWithSecond<DIM>, public GalerkinLpsIntegratorQ2<DIM>
{
  protected:

  public:

    std::string GetName() const {return "IntegratorLpsWithSecond";}
    void BasicInit() { GalerkinLpsIntegratorQ2<DIM>::BasicInit();}
}; 
}

/**********************************************************/

#endif
