/*----------------------------   stopwatch.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __stopwatch_H
#define __stopwatch_H
/*----------------------------   stopwatch.h     ---------------------------*/


/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2018 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include <chrono>
#include <ctime>
#include <cassert>
#include <string>
#include <map>
#include <iostream>
#include <fstream>

#include <iomanip>

/*----------------------------------------------------*/

namespace Gascoigne
{
  /**
   *
   * StopWatch that measure the 'real time'
   *
   **/
  class StopWatch 
  {
  protected:
    std::chrono::high_resolution_clock::time_point last_time;
    std::chrono::duration<double> sum_time;

    bool     running;
  
  public:
    StopWatch() : sum_time(std::chrono::duration<double>::zero()), running(false)
    {}

    virtual void  reset()
    {
      sum_time = std::chrono::duration<double>::zero();
      running = false;
    }
    virtual void  start()
    {
      assert(!running);
      last_time = std::chrono::high_resolution_clock::now();
    }
    virtual void stop()
    {
      std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
      running = false;
      sum_time += std::chrono::duration_cast<std::chrono::duration<double> >(now - last_time);
    }
    virtual double read() const
    {
      assert(!running);
      return sum_time.count();
    }
    virtual double read100() const
    {
      return static_cast<double>(static_cast<int>(read()*100))/100.0;
    }
  };

 
  /**
   *
   * StopWatch that measure the 'cpu time'
   *
   **/
  
  class CPUStopWatch : public StopWatch
  {
  protected:
    std::clock_t last_time;
    double sum_time;
  
  public:
    CPUStopWatch() : sum_time(0)
    {}

    void  reset()
    {
      sum_time = 0;
      running = false;
    }
    void  start()
    {
      assert(!running);
      last_time = std::clock();
    }
    void stop()
    {
      std::clock_t now = std::clock();
      running = false;
      sum_time += static_cast<double>(now - last_time)/CLOCKS_PER_SEC;
    }
    double read() const
    {
      assert(!running);
      return sum_time;
    }
    double read100() const
    {
      return static_cast<double>(static_cast<int>(read()*100))/100.0;
    }
  };

  /*----------------------------------------------------*/

  /**
   *
   * class for managing StopWatches
   *
   **/

  class Timer
  {
  protected:
    std::map<std::string, StopWatch> watches;
    std::map<std::string, CPUStopWatch> cpuwatches;

    std::map<std::string, size_t> counters;

  public:
    void  reset()
    {
      watches.clear();
      cpuwatches.clear();
      counters.clear();
    }
    void reset(const std::string& label)
    {
      auto it = watches.find(label);
      assert(it!=watches.end());
      watches.erase(it);
      
      auto cpuit = cpuwatches.find(label);
      assert(cpuit!=cpuwatches.end());
      cpuwatches.erase(cpuit);
      
      auto counter = counters.find(label);
      assert(counter!=counters.end());
      counters.erase(counter);
    }
    
    void  start(const std::string& label)
    {
      watches[label].start();
      cpuwatches[label].start();
    }

    void stop(const std::string& label)
    {
      assert(watches.find(label)!=watches.end());
      assert(cpuwatches.find(label)!=cpuwatches.end());
      watches[label].stop();
      cpuwatches[label].stop();
    }

    void count(const std::string& label){
      counters[label]++;
    }

    double read(const std::string& label) const
    {
      auto it = watches.find(label);
      assert(it!=watches.end());
      return it->second.read();
    }

    double read100(const std::string& label) const
    { 
      auto it = watches.find(label);
      assert(it!=watches.end());
      return it->second.read100();
    }

    double cpuread(const std::string& label) const
    {
      auto it = cpuwatches.find(label);
      assert(it!=cpuwatches.end());
      return it->second.read();
    }
    double cpuread100(const std::string& label) const
    { 
      auto it = cpuwatches.find(label);
      assert(it!=cpuwatches.end());
      return it->second.read100();
    }

    size_t counterread(const std::string& label) const {
      auto it = counters.find(label);
      assert(it!=counters.end());
      return it->second;
    }

    void print(const std::string& label) const
    {
      auto it = watches.find(label);
      assert(it!=watches.end());
      auto cit = cpuwatches.find(label);
      assert(cit!=cpuwatches.end());
      
      std::cout << std::setw(24) << std::left << label << std::setw(12) << std::right << it->second.read() << std::setw(12) << std::right << cit->second.read() << std::endl;
    }

    void print100(const std::string& label) const
    {
      auto it = watches.find(label);
      assert(it!=watches.end());
      auto cit = cpuwatches.find(label);
      assert(cit!=cpuwatches.end());

      std::cout << label << "\t" << it->second.read100() << "\t" << cit->second.read100() << std::endl;
    }
    
    void print100tofile(const std::string& filename) const
    {
      std::ofstream watch_logfile(filename);
      watch_logfile.precision(12);
      for (auto it : watches)
      {
        auto itt = watches.find(it.first);
        assert(itt!=watches.end());
        auto cit = cpuwatches.find(it.first);
        assert(cit!=cpuwatches.end());

        watch_logfile<< it.first << "\t" << itt->second.read100() << "\t" << cit->second.read100() << std::endl;
      }
      watch_logfile.close();
    }
    
    void print() const
    {
      std::cout << std::showpoint << std::fixed << std::setprecision(6);
      for (auto it : watches){
	      print(it.first);
      }
      for(auto it : counters){
        std::cerr << std::setw(24) << std::left << it.first << std::setw(12) << std::right << it.second << std::endl;
      }
      std::cout << std::endl << std::noshowpoint;
    }

    void print100() const
    {
      for (auto it : watches) {
	      print100(it.first);
      }
      std::cout << std::endl;
    }

  };
 
 
}



/*----------------------------   stopwatch.h     ---------------------------*/
/* end of #ifndef __stopwatch_H */
#endif
/*----------------------------   stopwatch.h     ---------------------------*/
