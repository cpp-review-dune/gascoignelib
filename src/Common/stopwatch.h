#ifndef STPWATCH_H
#define STPWATCH_H

#include  <ctime>
#include  "gostream.h"


/***********************************************************
 *****   ueberarbeitete und fehlerbereinigte Version   ***** 
 *****   Michael Schmich                  15.01.2002   *****
 ***********************************************************/

/*----------------------------------------------------*/

namespace Gascoigne
{
class Time
{
 protected:

  double sec;
  int min, hour, day;
  
 public:

  Time() : sec(0.), min(0), hour(0),day(0) {}

  void reset() {sec=0.; min=hour=day=0;}

  double GetSeconds () const {return sec;}
  int    GetMinutes () const {return min;}
  int    GetHours   () const {return hour;}
  int    GetDays    () const {return day;}

  double GetTotalSeconds () const;

  void add(double s);

  friend std::ostream& operator<<(std::ostream &s, const Time& A);
};

/*----------------------------------------------------*/

class StopWatch 
{
 protected:
  bool     running;
  clock_t  last_time;
  Time     T;
  
 public:
  StopWatch();

  const Time& GetTime() const {return T;}

  void  add(double s) { T.add(s);}
  void  reset();
  void  start();
  double stop() ;
  double read() const;
};
}

#endif
