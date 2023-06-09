/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#ifndef __BoundaryFunctional_h
#define __BoundaryFunctional_h

#include "functional.h"
#include "vertex.h"
#include <set>

/*-----------------------------------------*/

namespace Gascoigne {
class BoundaryFunctional : public virtual Functional
{
private:
protected:
public:
  BoundaryFunctional()
    : Functional()
  {}
  virtual ~BoundaryFunctional(){};

  virtual double J(const FemFunction& U,
                   const Vertex2d& v,
                   const Vertex2d& n,
                   int color) const
  {
    std::cerr << "\"BoundaryFunctional::J\" for 2d not written!" << std::endl;
    abort();
  }
  virtual double J(const FemFunction& U,
                   const Vertex3d& v,
                   const Vertex3d& n,
                   int color) const
  {
    std::cerr << "\"BoundaryFunctional::J\" for 3d not written!" << std::endl;
    abort();
  }
};
} // namespace Gascoigne

#endif
