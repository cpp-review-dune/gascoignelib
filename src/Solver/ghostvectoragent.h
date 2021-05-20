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

#ifndef __GhostVectorAgent_h
#define __GhostVectorAgent_h

#include "compvector.h"
#include "gascoigne.h"
#include "vectorinterface.h"
#include <string>

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVectorAgent
////
////
/////////////////////////////////////////////

class GhostVectorAgent : public std::map<Vector, GlobalVector*>
{
public:
  typedef std::map<Vector, GlobalVector*>::const_iterator const_iterator;
  typedef std::map<Vector, GlobalVector*>::iterator iterator;

  //
  ////  Con(De)structor
  //

  GhostVectorAgent();
  ~GhostVectorAgent();

  void Register(const Vector& mg);
  void Delete(Vector& mg);

  GlobalVector& operator()(const Vector& g);

  friend std::ostream& operator<<(std::ostream& os, const GhostVectorAgent& gva)
  {
    int i = 0, n = gva.size();
    os << "GhostVectorAgent: size=" << n << ", ";
    for (auto p = gva.begin(); p != gva.end(); p++, i++) {
      os << "VectorInterface(" << i << ")=('" << p->first.GetName() << "',"
         << p->second << ")";
      if (i < n - 1)
        os << ", ";
      else
        os << ". ";
    }
    return os;
  }
};
} // namespace Gascoigne

#endif
