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

  GhostAgent();
  ~GhostAgent();

  void Register(const std::string& mg);
  void Delete(std::string& mg);

  T& operator()(const std::string& g);
  const T& operator()(const std::string& g) const;

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
