/**
 *
 * Copyright (C) 2004 by the Gascoigne 3D authors
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

#ifndef __ilupermutate_H
#define __ilupermutate_H

#include "../Common/compvector.h"
#include "../DofHandler/dofhandlerbase.h"
#include "../LinAlg/columnstencil.h"

namespace Gascoigne {
class StreamDirection
{
  int dimension;
  int dx, dy, dz;
  const DofHandlerBase* M;
  const ColumnStencil* S;
  const GlobalVector& X;

  void Permutate(IndexVector& perm);

public:
  StreamDirection(const DofHandlerBase* m,
                  const StencilInterface* s,
                  const GlobalVector& x);

  void Permutate(IndexVector& perm, const IndexVector d);

  bool operator()(int i, int j) const;
  double est(int i, int j) const;
};

class VecDirection
{
  Vertex2d dir2d;
  Vertex3d dir3d;
  int dimension;
  const DofHandlerBase* M;

  void Permutate(IndexVector& perm);

public:
  VecDirection(const DofHandlerBase* m);

  void Permutate(IndexVector& perm, DoubleVector v);

  bool operator()(int i, int j) const;
};
} // namespace Gascoigne

#endif
