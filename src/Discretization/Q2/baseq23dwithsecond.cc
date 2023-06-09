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

#include "baseq23dwithsecond.h"

using namespace std;

namespace Gascoigne {

#define NDOF 27
#define NDOF1d 3
#define NDOF2d 9

/*-----------------------------------------*/

BaseQ23dWithSecond::BaseQ23dWithSecond()
  : BaseQ23d()
{
  dxx.reservesize(NDOF);
  dxy.reservesize(NDOF);
  dxz.reservesize(NDOF);
  dyy.reservesize(NDOF);
  dyz.reservesize(NDOF);
  dzz.reservesize(NDOF);
}

/*------------------------------------------------------------------*/

void
BaseQ23dWithSecond::point(const Vertex3d& s) const
{
  BaseQ23d::point(s);

  for (int i = 0; i < NDOF; i++) {
    int ix = (i % NDOF2d) % NDOF1d;
    int iy = (i % NDOF2d) / NDOF1d;
    int iz = i / NDOF2d;

    dxx[i] = psi_xx(ix, s.x()) * psi(iy, s.y()) * psi(iz, s.z());
    dxy[i] = psi_x(ix, s.x()) * psi_x(iy, s.y()) * psi(iz, s.z());
    dxz[i] = psi_x(ix, s.x()) * psi(iy, s.y()) * psi_x(iz, s.z());
    dyy[i] = psi(ix, s.x()) * psi_xx(iy, s.y()) * psi(iz, s.z());
    dyz[i] = psi(ix, s.x()) * psi_x(iy, s.y()) * psi_x(iz, s.z());
    dzz[i] = psi(ix, s.x()) * psi(iy, s.y()) * psi_xx(iz, s.z());
  }
}

/*------------------------------------------------------------------*/

#undef NDOF
#undef NDOF1d
#undef NDOF2d

} // namespace Gascoigne
