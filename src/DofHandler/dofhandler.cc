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

#include "dofhandler.h"

/*-----------------------------------------*/

namespace Gascoigne
{
template <>
std::string DofHandler<2>::GetName() const
{
  return "DofHandler 2d";
}
template <>
std::string DofHandler<3>::GetName() const
{
  return "DofHandler 3d";
}

template <>
IntVector DofHandler<2>::IndicesOfCell(int iq) const
{
  IntVector indices(4);

  int iq4 = iq * 4;

  indices[0] = nc[iq4];
  indices[1] = nc[iq4 + 1];
  indices[2] = nc[iq4 + 3];
  indices[3] = nc[iq4 + 2];

  return indices;
}
template <>
IntVector DofHandler<3>::IndicesOfCell(int iq) const
{
  IntVector indices(8);

  int offset = 8 * iq;

  indices[0] = nc[offset];
  indices[1] = nc[offset + 1];
  indices[2] = nc[offset + 3];
  indices[3] = nc[offset + 2];
  indices[4] = nc[offset + 4];
  indices[5] = nc[offset + 5];
  indices[6] = nc[offset + 7];
  indices[7] = nc[offset + 6];

  return indices;
}

template <>
const Vertex<2>& DofHandler<2>::vertex2d(int i) const
{
  return vertex(i);
}
template <>
const Vertex<3>& DofHandler<3>::vertex3d(int i) const
{
  return vertex(i);
}
// wrong dimension
template <>
const Vertex<3>& DofHandler<2>::vertex3d(int i) const
{
  std::cerr << "Error: no 3d-vertex in 2d!" << std::endl;
  abort();
}
template <>
const Vertex<2>& DofHandler<3>::vertex2d(int i) const
{
  std::cerr << "Error: no 2d-vertex in 3d!" << std::endl;
  abort();
}

}  // namespace Gascoigne

template class Gascoigne::DofHandler<2>;
template class Gascoigne::DofHandler<3>;
