#ifndef  __Timer_h
#define  __Timer_h

#include  <map>
#include  <iostream>

#include  "stopwatch.h"
#include  "gascoigne.h"
#include  <string>
#include  <cassert>

/*-----------------------------------------*/

namespace Gascoigne
{
class Timer
{
protected:

  typedef std::map<std::string,StopWatch>::const_iterator  const_iterator;
  typedef std::map<std::string,StopWatch>::iterator        iterator;

  std::map<std::string,StopWatch>  mstpw;

public:
 
  Timer();
  Timer(const Timer& T)
    {
      const_iterator p = T.begin();
      for(;p!=T.end();p++)
	{
	  mstpw.insert(*p);
	}
    }

  void init(const Timer& T)
    {
      assert(0);
      for(const_iterator p = T.begin(); p!= T.end();p++)      add(p->first);
    }

  void add(const Timer& T)
    {
      for(iterator p = mstpw.begin(); p!= mstpw.end();p++)
	{
	  double s = T.read(p->first);
	  p->second.add(s);
	}
    }
  
  void add(const std::string& s) { mstpw.insert(std::make_pair(s,StopWatch() ));}

  const StopWatch& Get(const std::string& s) const
    {
      const_iterator p =  mstpw.find(s);
      if(p==mstpw.end()) {std::cerr << "no such StopWatch: "<<s<<std::endl; abort();}
      return p->second;
    }

  const_iterator begin() const {return mstpw.begin();}
  const_iterator end  () const {return mstpw.end  ();}
  void start(const std::string& s)
    {
      iterator p = mstpw.find(s);
      if(p==mstpw.end())
	{
	  std::cerr << "Timer::start()\n";
	  std::cerr << "no such StopWatch: " << s << std::endl; abort();
	}
      p->second.start();
    }
  void stop(const std::string& s)
    {
      iterator p = mstpw.find(s);
      if(p==mstpw.end())
	{
	  std::cerr << "Timer::stop()\n";
	  std::cerr << "no such StopWatch: " << s << std::endl; abort();
	}
      p->second.stop();
    }
  double read(const std::string& s) const
    {
      const_iterator p = mstpw.find(s);
      if(p==mstpw.end())
	{
	  std::cerr << "Timer::read()\n";
	  std::cerr << "no such StopWatch: " << s << std::endl; abort();
	}
      return p->second.read();
    }
  StopWatch total() const 
    {
      StopWatch s;
      const_iterator p = mstpw.begin();
      while(p!=mstpw.end()) {s.add(p++->second.read());}
      return s;
    }
  friend std::ostream& operator<<(std::ostream& os, const Timer& T);
};
}


#endif
