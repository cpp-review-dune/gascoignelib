/**
 *
 * Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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

#ifndef __index_h
#define __index_h

#include <map>
#include <set>

#include "../Common/stlio.h"
#include "../Interface/gascoigne.h"

/*------------------------------------------*/

namespace Gascoigne {
class Index
{
protected:
  IndexVector vl2g, el2g, hl2g, ql2g;
  IndexMap vg2l, eg2l, hg2l, qg2l;

public:
protected:
  ////////////////
  // Local To Global
  ////////////////

  const IndexVector& Vertexl2g() const { return vl2g; }
  const IndexVector& Edgel2g() const { return el2g; }
  IndexVector& Edgel2g() { return el2g; }

  int Edgel2g(int i) const { return el2g[i]; }

public:
  ////////////////
  // Constructors
  ////////////////

  Index();
  Index(const Index& I);
  Index& operator=(const Index& I);

  int Vertexl2g(int i) const { return vl2g[i]; }

  IndexVector& Quadl2g() { return ql2g; }
  IndexVector& Hexl2g() { return hl2g; }
  const IndexVector& Hexl2g() const { return hl2g; }
  const IndexVector& Quadl2g() const { return ql2g; }
  int Quadl2g(size_t i) const
  {
    assert(i < ql2g.size());
    return ql2g[i];
  }
  int Hexl2g(size_t i) const
  {
    assert(i < hl2g.size());
    return hl2g[i];
  }
  const IndexMap& Quadg2l() const { return qg2l; }
  const IndexMap& Hexg2l() const { return hg2l; }

  ////////////////
  // Sizes
  ////////////////

  size_t nnodes() const { return VertexSize(); }
  size_t VertexSize() const { return vl2g.size(); }
  size_t VertexGSize() const { return vg2l.size(); }
  size_t EdgeSize() const { return el2g.size(); }
  size_t EdgeGSize() const { return eg2l.size(); }
  size_t HexSize() const { return hl2g.size(); }
  size_t HexGSize() const { return hg2l.size(); }
  size_t QuadSize() const { return ql2g.size(); }
  size_t QuadGSize() const { return qg2l.size(); }

  ////////////////
  // Global To Local
  ////////////////

  const IndexMap& Vertexg2l() const { return vg2l; }
  const IndexMap& Edgeg2l() const { return eg2l; }
  IndexMap& Edgeg2l() { return eg2l; }

  int Quadg2l(int i) const
  {
    const auto ip = qg2l.find(i);
    if (ip == qg2l.end()) {
      std::cerr << "Index:: Quadg2l" << std::endl;
      std::cerr << "there is no "
                << " " << i << std::endl;
      abort();
    }
    return ip->second;
  }
  int Hexg2l(int i) const
  {
    const auto ip = hg2l.find(i);
    if (ip == hg2l.end()) {
      std::cerr << "Index:: Hexg2l" << std::endl;
      std::cerr << "there is no "
                << " " << i << std::endl;
      abort();
    }
    return ip->second;
  }
  int Vertexg2l(int i) const
  {
    const auto ip = vg2l.find(i);
    if (ip == vg2l.end()) {
      std::cerr << "Index:: Vertexg2l" << std::endl;
      std::cerr << "there is no "
                << " " << i << std::endl;
      abort();
    }
    return ip->second;
  }
  int Edgeg2l(int i) const
  {
    const auto ip = eg2l.find(i);
    if (ip == eg2l.end()) {
      std::cerr << "Index:: Edgeg2l" << std::endl;
      std::cerr << "there is no "
                << " " << i << std::endl;
      abort();
    }
    return ip->second;
  }

  ////////////////////////
  // Global To Local Check
  ////////////////////////

  // gibt -2 zurueck, falls globaler vertex nicht in levelmesh

  int Vertexg2lCheck(int i) const
  {
    const auto ip = vg2l.find(i);
    if (ip == vg2l.end())
      return -2;
    return ip->second;
  }
  int Edgeg2lCheck(int i) const
  {
    const auto ip = eg2l.find(i);
    if (ip == eg2l.end())
      return -2;
    return ip->second;
  }
  int Quadg2lCheck(int i) const
  {
    const auto ip = qg2l.find(i);
    if (ip == qg2l.end())
      return -2;
    return ip->second;
  }
  int Hexg2lCheck(int i) const
  {
    const auto ip = hg2l.find(i);
    if (ip == hg2l.end())
      return -2;
    return ip->second;
  }

  friend std::ostream& operator<<(std::ostream& os, const Index& I);

  void InitNodes(const IndexSet& nodes);
  void InitEdges(const IndexSet& edges);
  void InitQuads();
  void InitHexs();
};
} // namespace Gascoigne

#endif
