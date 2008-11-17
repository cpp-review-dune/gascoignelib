#include "stopwatch.h"

//#include  <strstream>
#include  <stdlib.h>
#include  <stdio.h>
#include  <sys/times.h>
#include  <sys/types.h>
#include  <sys/wait.h>
#include  <unistd.h>

using namespace std;

/***********************************************************
 *****   ueberarbeitete und fehlerbereinigte Version   ***** 
 *****   Michael Schmich                  15.01.2002   *****
 ***********************************************************/

/*----------------------------------------------------*/

namespace Gascoigne
{

double Time::GetTotalSeconds () const 
{
  return sec+60.*min+3600.*hour+24.*3600.*day;
}

/*----------------------------------------------------*/

void Time::add(double s)
{
  sec += s;
  int su = static_cast<int>(sec)/60;

  sec -= 60.*su;
  min += su;

  hour += min / 60;
  min   = min % 60;

  day += hour / 24;  
  hour = hour % 24;

}

/*----------------------------------------------------*/

ostream& operator<<(ostream &s, const Time& A)
{
  if(A.GetDays()>0)
    {
	s.precision(1);
	s << A.GetDays() << "-";
    }
  s.precision(2);
  if(A.GetHours()>0)
    {
	s << A.GetHours() << "-";
    }
  s << A.GetMinutes() << ":" << A.GetSeconds();
  return s;
}

/*----------------------------------------------------*/

StopWatch::StopWatch() : running(0), last_time(0), T() {}

/*----------------------------------------------------*/

void StopWatch::reset() 
{ 
  running = 0; last_time = 0; T.reset(); 
}

/*----------------------------------------------------*/

void StopWatch::start() 
{ 
  if (!running) { last_time = clock(); running = 1;}
}

/*----------------------------------------------------*/

double StopWatch::stop()  
{ 
  if (running) 
    {
	double s = static_cast<double>(clock() - last_time) / static_cast<double>(CLOCKS_PER_SEC);
	add(s);
	running = 0;
    }
  return T.GetTotalSeconds(); 
}


/*----------------------------------------------------*/

double StopWatch::read() const  
{
  if (running) return -1;
  return T.GetTotalSeconds(); 
} 


double StopWatch::read100() const  
{
  if (running) return -1;
  double X = static_cast<double> (0.01 * static_cast<int> (100 * read()));
  return X;

}

}
