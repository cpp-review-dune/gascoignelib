#ifndef  __Application_h
#define  __Application_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Application

////
////
/////////////////////////////////////////////

#include  "gascoigne.h"


/*-------------------------------------------------------*/

namespace Gascoigne
{
  class Application
  {
    private:
      mutable double _dt, _time;

    protected:
      double GetTime() const {
        return _time;
      }
      double GetTimeStep() const {
        return _dt;
      }

    public:
      //
      ////  Con(De)structor 
      //

      Application() : _dt(0.),_time(0.) {}
      virtual ~Application() {}

      virtual std::string GetName() const=0;

      virtual void SetTime(double time, double dt) const {
        _time = time; 
        _dt = dt;
      }

      virtual void SetFemData(FemData& q) const {}
      virtual void SetParameterData(LocalParameterData& q) const {}
  };
}

#endif

