/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#ifndef __cell_h
#define __cell_h

#include "gascoigne.h"
#include "stlio.h"
#include "vertex.h"
#include <array>
#include <string>

/*-----------------------------------------------------------*/

// N = no of nodes
// E = no of edges/faces

/*-----------------------------------------------------------*/

namespace Gascoigne {
template<size_t N, int E>
class Cell : public std::array<IndexType, N> /* das sind die vertex-no. */
{
protected:
  /* Data */

  IndexType qlevel, qfather, mat, mat_Vanka;
  std::array<Vertex3d, 3> bas_Vanka;
  IndexVector qchilds;
  std::array<IndexType, E> qedges; /* edge numbers */

public:
  /* Constructors */

  Cell()
    : std::array<IndexType, N>()
    , qlevel(0)
    , qfather(-1)
    , mat(0)
    , mat_Vanka(0)
  {
    qedges.fill(-1);
  }

  Cell(const Cell& c)
    : std::array<IndexType, N>(c)
    , qlevel(c.level())
    , qfather(c.father())
    , mat(c.material())
    , mat_Vanka(c.material_Vanka())
    , bas_Vanka(c.basis_Vanka())
    , qchilds(c.childs())
    , qedges(c.edges())
  {}

  Cell(IndexType l, IndexType f)
    : std::array<IndexType, N>()
    , qlevel(l)
    , qfather(f)
    , mat(0)
    , mat_Vanka(0)
    , qchilds(0)
    , qedges()
  {
    this->fill(-17);
  }

  /* Operators */

  Cell<N, E>& operator=(const Cell<N, E>& c)
  {
    vertex() = c.vertex();
    qlevel = c.level();
    qfather = c.father();
    qchilds.memory(c.nchilds());
    qchilds = c.childs();
    qedges = c.edges();
    return *this;
  }
  bool operator==(const Cell<N, E>& c) const
  {
    if (vertex() == c.vertex())
      return 1;
    return 0;
  }

  /* Zugriff */

  IndexType level() const { return qlevel; }
  IndexType& level() { return qlevel; }
  IndexType father() const { return qfather; }
  IndexType& father() { return qfather; }
  bool sleep() const { return qchilds.size() != 0; }
  IndexType nchilds() const { return qchilds.size(); }
  IndexType nvertexs() const { return N; }
  IndexType child(IndexType i) const { return qchilds[i]; }
  IndexType& child(IndexType i) { return qchilds[i]; }
  IndexType vertex(IndexType i) const { return (*this)[i]; }
  IndexType& vertex(IndexType i) { return (*this)[i]; }
  IndexType edge(IndexType i) const { return qedges[i]; }
  IndexType& edge(IndexType i) { return qedges[i]; }

  const IndexVector& childs() const { return qchilds; }
  IndexVector& childs() { return qchilds; }
  const std::array<IndexType, N>& vertex() const { return (*this); }
  std::array<IndexType, N>& vertex() { return (*this); }
  const std::array<IndexType, E>& edges() const { return qedges; }
  std::array<IndexType, E>& edges() { return qedges; }

  IndexType material() const { return mat; }
  IndexType& material() { return mat; }

  IndexType material_Vanka() const { return mat_Vanka; }
  IndexType& material_Vanka() { return mat_Vanka; }

  std::array<Vertex3d, 3> basis_Vanka() const { return bas_Vanka; }
  std::array<Vertex3d, 3>& basis_Vanka() { return bas_Vanka; }

  /* Functions */

  template<int M>
  void vertex_loc2glob(std::array<IndexType, M>& ig,
                       const std::array<IndexType, M>& il) const
  {
    typename std::array<IndexType, M>::iterator gp = ig.begin();
    typename std::array<IndexType, M>::const_iterator lp = il.begin();
    while (lp != il.end())
      *gp++ = (*this)[*lp++];
  }

  IndexType global2local(IndexType gi) const;

  void BinWrite(std::ostream& s) const
  {
    ArrayBinWrite(s, vertex());
    int sizeInt = sizeof(int);
    s.write(reinterpret_cast<const char*>(&qlevel), sizeInt);
    s.write(reinterpret_cast<const char*>(&qfather), sizeInt);
    IndexType nc = nchilds();
    s.write(reinterpret_cast<const char*>(&nc), sizeInt);
    for (IndexType i = 0; i < nchilds(); i++) {
      s.write(reinterpret_cast<const char*>(&qchilds[i]), sizeInt);
    }
    ArrayBinWrite(s, edges());
  }

  void BinRead(std::istream& s)
  {
    ArrayBinRead(s, vertex());
    IndexType sizeInt = sizeof(int);
    s.read(reinterpret_cast<char*>(&qlevel), sizeInt);
    s.read(reinterpret_cast<char*>(&qfather), sizeInt);
    IndexType nc;
    s.read(reinterpret_cast<char*>(&nc), sizeInt);
    childs().resize(nc);
    for (IndexType i = 0; i < nchilds(); i++) {
      s.read(reinterpret_cast<char*>(&qchilds[i]), sizeInt);
    }
    ArrayBinRead(s, edges());
  }

  friend std::ostream& operator<<(std::ostream& s, const Cell& A)
  {
    s << A.vertex() << " ";
    s << A.level() << " ";
    s << A.father() << " ";
    s << A.material() << " @ ";
    s << A.nchilds() << " " << A.childs();
    s << " : " << A.edges();
    s << std::endl;

    return s;
  }

  friend std::istream& operator>>(std::istream& s, Cell& A)
  {
    std::string symbol;
    IndexType n;
    s >> A.vertex();
    s >> A.level();
    s >> A.father();

    // old version without material?
    bool oldversion = false;
    s >> symbol;
    try {
      A.material() = std::stoi(symbol);
    } catch (const std::exception& e) {
      A.material() = 0;
      oldversion = true;
    }
    if (!oldversion) // new version: mat was correct,
      s >> symbol;   // now read symbol
    //

    if (symbol != "@") {
      std::cout << "ERROR: Cell::operator>>" << std::endl;
      exit(1);
    }

    s >> n;
    A.childs().resize(n);
    s >> A.childs();
    s >> symbol;
    if (symbol != ":") {
      std::cout << "ERROR: Cell::operator>>" << std::endl;
      exit(1);
    }
    s >> A.edges();

    return s;
  }
};

template<size_t N, int E>
inline IndexType
Cell<N, E>::global2local(IndexType gi) const
{
  for (IndexType i = 0; i < N; i++) {
    if (vertex(i) == gi)
      return i;
  }
  return -1;
}
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
