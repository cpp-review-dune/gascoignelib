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


#ifndef __Pi_h
#define __Pi_h

#include "compvector.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include <map>
#include <array>

namespace Gascoigne
{
  /*-----------------------------------------*/

  class Pi
  {
  protected:
    std::map<int, std::array<int, 2>> edge;
    std::map<int, std::array<int, 4>> face;
    std::map<int, std::array<int, 8>> cell;

    void Init2d(const GascoigneMesh2d *MP);
    void Init3d(const GascoigneMesh3d *MP);

  public:
    Pi();

    void Init(const GascoigneMesh *MP);

    void vmult(CompVector<double> &y,
               const CompVector<double> &x,
               double s = 1.) const;
  };
}

#endif
