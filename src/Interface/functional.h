#ifndef  __Functional_h
#define  __Functional_h

#include  "application.h"


/*-----------------------------------------*/

namespace Gascoigne
{

  //////////////////////////////////////////////
  //
  /// Interface class for Functional
  //
  //////////////////////////////////////////////

  class Functional : public virtual Application
  {
    private:
      
    protected:
      double  exact;
      bool    exactisknown;

    public:
      Functional() : Application(), exactisknown(0), exact(0.) {}
      ~Functional() {}
      Functional(const Functional& F) : Application(F) {
        exact = F.ExactValue();
      } 

      double  ExactValue() const { 
        return exact;
      }
      double& ExactValue() {
        exactisknown = 1; 
        return exact;
      }
      bool ExactValueIsKnown() const 
      { 
        return exactisknown; 
      }
  };
}

#endif
