/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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

#ifndef __localhierarchicalmesh3d_h
#define __localhierarchicalmesh3d_h

#include "hierarchicalmesh.h"
#include "hierarchicalmesh3d.h"

/*---------------------------------------------------*/

namespace Gascoigne {
class LocalHierarchicalMesh3d : public HierarchicalMesh3d {
protected:
  /*  typedef  */

  void inner_vertex_newton3d(const IntVector &, const IntSet &, const IntSet &);

public:
  void write_gup(const std::string &, const GlobalVector &u) const;
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
