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

class Application
{
private:

  mutable double _dt, _time;

protected:

  double GetTime() const {return _time;}
  double GetTimeStep() const {return _dt;}

public:

//
////  Con(De)structor 
//

  Application() : _dt(0.),_time(0.) {}
  virtual ~Application() {}
  Application(const Application&) {}

  virtual std::string GetName() const=0;

  virtual void SetTime(double time) const { _time = time;}
  virtual void SetTime(double time, double dt) const { _time = time; _dt = dt;}

  virtual void SetFemData(Gascoigne::FemData& q) const {}
  virtual void SetParameterData(Gascoigne::LocalData& q) const {}

  virtual int GetNcomp() const {assert(0);}

  int size()const{std::cerr<< "Application: never use \"size()\"\n"; assert(0);}
  int ncomp()const{std::cerr<< "Application: never use \"ncomp()\"\n"; assert(0);}
};


#endif


