#ifndef  __TimeInfo_h
#define  __TimeInfo_h

#include  <string>

/*-----------------------------------------*/

class TimeInfo
{
protected:

  double _deltat, _time, _theta;
  double _tbegin, _tend;
  int    _iter, _neuler;
  double _ftscale[3], _fttheta[3];
  std::string _scheme;
  double _Theta;

  int ftstep() const;

public:

  TimeInfo();

  void ReInit();
  void ReInit(double det);
  void BasicInit();

  double  dt   () const;
  double theta() const;
  void iteration(int i);
  void ReInit(double tb, double det, double te, const std::string& sch, int ne, double t);
  void scale_timestep(double s) { _deltat *= s;}
  void stepback() { _time -= _deltat;}
  double time() const { return _time;}

  double time_begin() const { return _tbegin;}
  double time_end()   const { return _tend;}

  void ReInitBackward(int niter, double endtime);
  void iteration_backward(int i);
};

#endif
