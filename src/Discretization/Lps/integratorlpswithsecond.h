#ifndef __IntegratorLpsWithSecond_h
#define __IntegratorLpsWithSecond_h

#include "integratorwithsecond.h"
#include "galerkinlpsintegratorq2.h"

/**********************************************************/

template<int DIM>
class IntegratorLpsWithSecond : public IntegratorWithSecond<DIM>, public Gascoigne::GalerkinLpsIntegratorQ2<DIM>
{
  protected:

  public:

    std::string GetName() const {return "IntegratorLpsWithSecond";}
}; 

/**********************************************************/

#endif
