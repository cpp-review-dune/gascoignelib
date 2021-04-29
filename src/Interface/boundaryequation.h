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

#ifndef __BoundaryEquation_h
#define __BoundaryEquation_h

#include "application.h"
#include "entrymatrix.h"
#include "compvector.h"
#include "vertex.h"

namespace Gascoigne {

/**********************************************************/

class BoundaryEquation : public virtual Application {
private:
protected:
public:
  BoundaryEquation() : Application() {}
  ~BoundaryEquation() {}

  /**
     clones a BoundaryEquation. Usually it is simply to return a new instance of
     the same object passing the required variables to the constructor. It takes
     the role of a copy constructor and the cloning of classes is required for
     multithreading.
  */
  virtual BoundaryEquation *createNew() const {
    std::cerr << "\"BoundaryEquation::createNew\" not written!" << std::endl;
    abort();
  }

  virtual int GetNcomp() const = 0;

  virtual void Form(VectorIterator b, const FemFunction &U,
                    const TestFunction &N, int col) const = 0;
  virtual void Matrix(EntryMatrix &E, const FemFunction &U,
                      const TestFunction &M, const TestFunction &N,
                      int col) const = 0;

  virtual void pointboundary(double h, const FemFunction &U, const Vertex2d &v,
                             const Vertex2d &n) const {}
  virtual void pointboundary(double h, const FemFunction &U, const Vertex3d &v,
                             const Vertex3d &n) const {}
  virtual void pointmatrixboundary(double h, const FemFunction &U,
                                   const Vertex2d &v, const Vertex2d &n) const {
    pointboundary(h, U, v, n);
  }
  virtual void pointmatrixboundary(double h, const FemFunction &U,
                                   const Vertex3d &v, const Vertex3d &n) const {
    pointboundary(h, U, v, n);
  }
};

/**********************************************************/

} // namespace Gascoigne

#endif
