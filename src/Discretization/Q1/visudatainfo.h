/**
 *
 * Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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

#ifndef __visudatainfo_h
#define __visudatainfo_h

#include <array>

#include "gascoigne.h"
#include "visudata.h"

#include <map>
#include <string>

/*-------------------------------------------------------------------------*/

namespace Gascoigne {
class VisuDataInfo
{
protected:
  std::map<std::string, IndexType> scalars;
  std::map<std::string, std::array<int, 3>> vectors;
  std::map<std::string, IndexType> scalar_order;
  std::map<std::string, IndexType> vector_order;

public:
  typedef std::map<std::string, IndexType>::const_iterator siterator;
  typedef std::map<std::string, std::array<int, 3>>::const_iterator viterator;

  VisuDataInfo() {}
  VisuDataInfo(IndexType ncomp) { AddScalars(ncomp); }
  VisuDataInfo(const VisuData& V, std::string def = "U");
  VisuDataInfo(const VisuDataInfo& V)
    : scalars(V.Scalars())
    , vectors(V.Vectors())
  {}
  VisuDataInfo& operator=(const VisuDataInfo& V);

  bool operator!=(const VisuDataInfo& V) const;

  void Clear()
  {
    scalar_order.clear();
    vector_order.clear();
    scalars.clear();
    vectors.clear();
  }

  siterator GetSIterator(IndexType i)
  {
    for (siterator p = sbegin(); p != send(); p++) {
      std::string s = p->first;
      if (scalar_order[s] == i)
        return p;
    }
    abort();
  }
  viterator GetVIterator(IndexType i)
  {
    for (viterator p = vbegin(); p != vend(); p++) {
      std::string s = p->first;
      if (vector_order[s] == i)
        return p;
    }
    abort();
  }

  void AddScalar(IndexType index, const std::string& name, IndexType i)
  {
    scalar_order[name] = index;
    scalars[name] = i;
  }
  void AddVector(IndexType index,
                 const std::string& name,
                 const std::array<int, 3>& i)
  {
    vector_order[name] = index;
    vectors[name] = i;
  }

  void AddScalars(IndexType ncomp, std::string def = "U");

  IndexType nscalars() const { return scalars.size(); }
  IndexType nvectors() const { return vectors.size(); }

  const std::map<std::string, IndexType>& Scalars() const { return scalars; }
  const std::map<std::string, std::array<int, 3>>& Vectors() const
  {
    return vectors;
  }

  siterator sbegin() const { return scalars.begin(); }
  siterator send() const { return scalars.end(); }
  viterator vbegin() const { return vectors.begin(); }
  viterator vend() const { return vectors.end(); }
};
} // namespace Gascoigne

#endif
