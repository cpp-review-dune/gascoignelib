/**
 *
 * Copyright (C) 2004, 2005, 2020 by the Gascoigne 3D authors
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

#ifndef __GhostAgent_h
#define __GhostAgent_h

#include "gascoigne.h"
#include "matrixinterface.h"
#include <string>

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVectorAgent
////
////
/////////////////////////////////////////////

template<typename T>
class GhostAgent : public std::map<std::string, T*>
{
public:
  // typedef std::map<std::string, T*>::const_iterator const_iterator;
  // typedef std::map<std::string, T*>::iterator iterator;

  //
  ////  Con(De)structor
  //

  GhostAgent(){};
  ~GhostAgent()
  {
    for (auto p = std::map<std::string, T*>::begin();
         p != std::map<std::string, T*>::end();
         p++) {
      if (p->second) {
        delete p->second;
        p->second = NULL;
      }
    }
  };

  void Delete(std::string& mg)
  {
    auto p = std::map<std::string, T*>::find(mg);
    if (p != std::map<std::string, T*>::end()) {
      delete p->second;
      std::map<std::string, T*>::erase(p);
    }
  };

  T& operator()(const std::string& g)
  {
    auto p = std::map<std::string, T*>::find(g);
    if (p == std::map<std::string, T*>::end()) {
      throw std::runtime_error("Element not in GhostAgent: " + g);
    }
    T* vp = p->second;
    if (vp == NULL) {
      throw std::runtime_error("Element in GhostAgent not initialized: " + g);
    }
    return *vp;
  };

  const T& operator()(const std::string& g) const
  {
    auto p = std::map<std::string, T*>::find(g);
    if (p == std::map<std::string, T*>::end()) {
      throw std::runtime_error("Element not in GhostAgent: " + g);
    }
    T* vp = p->second;
    if (vp == NULL) {
      throw std::runtime_error("Element in GhostAgent not initialized: " + g);
    }
    return *vp;
  };

  friend std::ostream& operator<<(std::ostream& os, const GhostAgent& gva)
  {
    int i = 0, n = gva.size();
    os << "GhostAgent: size=" << n << ", ";
    for (auto p = gva.begin(); p != gva.end(); p++, i++) {
      os << "Vector(" << i << ")=('" << p->first << "'," << p->second << ")";
      if (i < n - 1)
        os << ", ";
      else
        os << ". ";
    }
    return os;
  }
};

using GhostVectorAgent = GhostAgent<GlobalVector>;
using MatrixAgent = GhostAgent<MatrixInterface>;

} // namespace Gascoigne

#endif
