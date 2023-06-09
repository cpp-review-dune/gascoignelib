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

#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <iomanip>

/*----------------------------------------------------*/

namespace Gascoigne {
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

  std::clock_t last_cpu_time;
  double sum_cpu_time;

  size_t activations;

  bool running;

public:
  StopWatch()
    : sum_time(std::chrono::duration<double>::zero())
    , sum_cpu_time(0)
    , activations(0)
    , running(false)
  {}

  virtual void reset()
  {
    sum_time = std::chrono::duration<double>::zero();
    sum_cpu_time = 0;
    running = false;
  }

  virtual void start()
  {
    if (!running) {
      activations++;
      last_time = std::chrono::high_resolution_clock::now();
      last_cpu_time = std::clock();
      running = true;
    }
  }

  virtual void stop()
  {
    std::chrono::high_resolution_clock::time_point now =
      std::chrono::high_resolution_clock::now();
    sum_time += std::chrono::duration_cast<std::chrono::duration<double>>(
      now - last_time);
    last_time = std::chrono::high_resolution_clock::now();

    std::clock_t now_cpu = std::clock();
    sum_cpu_time +=
      static_cast<double>(now_cpu - last_cpu_time) / CLOCKS_PER_SEC;
    last_cpu_time = std::clock();

    running = false;
  }

  virtual double read() const
  {
    assert(!running);
    return sum_time.count();
  }

  double readPerActivations() const
  {
    assert(!running);
    return sum_time.count() / activations;
  }

  virtual double readCPU() const
  {
    assert(!running);
    return sum_cpu_time;
  }

  virtual double read100() const
  {
    return static_cast<double>(static_cast<int>(read() * 100)) / 100.0;
  }

  double readCPU100() const
  {
    return static_cast<double>(static_cast<int>(read() * 100)) / 100.0;
  }

  size_t readActivations() const { return activations; }

  bool isRunning() const { return running; }
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
  std::map<std::string, size_t> counters;

public:
  const std::map<std::string, StopWatch>& GetWatches() const { return watches; }
  const std::map<std::string, size_t>& GetCounters() const { return counters; }

  void reset()
  {
    watches.clear();
    counters.clear();
  }
  void reset(const std::string& label)
  {
    auto it = watches.find(label);
    assert(it != watches.end());
    it->second.reset();

    auto counter = counters.find(label);
    assert(counter != counters.end());
    counters.erase(counter);
  }

  void start(const std::string& label) { watches[label].start(); }

  void stop(const std::string& label)
  {
    if (watches.find(label) == watches.end()) {
      std::cerr << "StopWatch: did not find running StopWatch with label "
                << label << std::endl;
      abort();
    }
    watches[label].stop();
  }

  void count(const std::string& label) { counters[label]++; }

  double read(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());
    return it->second.read();
  }

  double read100(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());
    return it->second.read100();
  }

  double cpuread(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());
    return it->second.readCPU();
  }
  double cpuread100(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());
    return it->second.readCPU100();
  }

  size_t counterread(const std::string& label) const
  {
    auto it = counters.find(label);
    assert(it != counters.end());
    return it->second;
  }

  void print(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());

    std::cout << std::setw(24) << std::left << label << std::setw(12)
              << std::right << it->second.read() << std::setw(12) << std::right
              << it->second.readCPU() << std::setw(12) << std::right
              << it->second.readPerActivations() << std::setw(12) << std::right
              << it->second.readActivations() << std::endl;
  }

  void print100(const std::string& label) const
  {
    auto it = watches.find(label);
    assert(it != watches.end());

    std::cout << label << "\t" << it->second.read100() << "\t"
              << it->second.readCPU100() << std::endl;
  }

  void print100tofile(const std::string& filename) const
  {
    std::ofstream watch_logfile(filename);
    watch_logfile.precision(12);
    for (const auto& it : watches) {
      watch_logfile << it.first << "\t" << it.second.read100() << "\t"
                    << it.second.readCPU100() << std::endl;
    }
    watch_logfile.close();
  }

  void print() const
  {
    std::cout << "\n*** Timer ***\n";
    std::cout << std::showpoint << std::fixed << std::setprecision(6);
    for (auto it : watches) {
      print(it.first);
    }
    for (auto it : counters) {
      std::cerr << std::setw(24) << std::left << it.first << std::setw(12)
                << std::right << it.second << std::endl;
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

} // namespace Gascoigne

/*----------------------------   stopwatch.h     ---------------------------*/
/* end of #ifndef __stopwatch_H */
#endif
/*----------------------------   stopwatch.h     ---------------------------*/
