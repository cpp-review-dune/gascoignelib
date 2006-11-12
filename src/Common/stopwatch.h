#ifndef STPWATCH_H
#define STPWATCH_H

#include  <ctime>
#include  "gostream.h"
#include  <unistd.h>
#include <sys/resource.h>
#include <time.h>

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

/*----------------------------------------------------*/

#ifdef __WITH_RSTStopWatch__
	  // measures the walltime
class RTStopWatch 
{
 protected:
  bool     running;
  timespec t1,t2;
  timespec my_time;
  
 public:
  RTStopWatch() 
    {
      reset();
    }
  

  void  reset()
    {
      running         = 0;
      my_time.tv_sec  = 0;
      my_time.tv_nsec = 0;
    }
  
  void  start()
    {
      if (!running) clock_gettime(CLOCK_REALTIME,&t1);
      running = 1;
    }
  
  double stop()
    {
      if (running)
	{
	  clock_gettime(CLOCK_REALTIME,&t2);
	  my_time.tv_sec  += t2.tv_sec-t1.tv_sec;
	  my_time.tv_nsec += t2.tv_nsec-t1.tv_nsec;

	  while (long(my_time.tv_nsec)>long(1000000000))
	    {
	      my_time.tv_nsec -= long(1000000000);
	      my_time.tv_sec  += 1;
	    }
	  while (long(my_time.tv_nsec)<long(0))
	    {
	      my_time.tv_nsec += long(1000000000);
	      my_time.tv_sec  -= 1;
	    }
	}
      running = 0;
      return double( my_time.tv_sec)+1.e-9 * double (my_time.tv_nsec);
    }
  
  double read() const
    {
      if (running) return -1;
      return double (my_time.tv_sec)+1.e-9 * double (my_time.tv_nsec);
    }
  
};

#endif /*__WITH_RSTStopWatch__*/


}

#endif
