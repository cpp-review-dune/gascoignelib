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

#include "visudatacompvector.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {
VisuDataCompVector::VisuDataCompVector()
  : VisuData()
  , _v(NULL)
{}

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector(const GlobalVector& v)
  : VisuData()
  , _v(NULL)
{
  SetGlobalVector(&v);
}

/*---------------------------------------------------*/

void
VisuDataCompVector::SetGlobalVector(const GlobalVector* v)
{
  assert(v);
  _v = v;
}

/*---------------------------------------------------*/

int
VisuDataCompVector::visucomp() const
{
  return _v->ncomp();
}

/*---------------------------------------------------*/

int
VisuDataCompVector::visun() const
{
  return _v->n();
}

/*---------------------------------------------------*/

double
VisuDataCompVector::visudata(int i, int c) const
{
  assert(_v);
  assert(i < _v->n());
  assert(c < _v->ncomp());
  return (*_v)(i, c);
}
} // namespace Gascoigne
