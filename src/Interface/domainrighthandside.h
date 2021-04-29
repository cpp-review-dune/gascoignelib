/**
 *
 * Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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

#ifndef __DomainRightHandSide_h
#define __DomainRightHandSide_h

#include "application.h"
#include "vertex.h"
#include "compvector.h"

namespace Gascoigne {

/*-------------------------------------------------------*/

class DomainRightHandSide : public virtual Application {
private:
protected:
public:
  DomainRightHandSide() {}
  ~DomainRightHandSide() {}

  /**
     clones a DRHS. Usually it is simply to return a new instance of the same
     object passing the required variables to the constructor. It takes the role
     of a copy constructor and the cloning of classes is required for
     multithreading.
  */
  virtual DomainRightHandSide *createNew() const {
    std::cerr << "\"DomainRightHandSide::createNew\" not written!" << std::endl;
    abort();
  }

  virtual int GetNcomp() const = 0;

  virtual double operator()(int c, const Vertex2d &v) const {
    std::cerr << "\"DomainRightHandSide::operator()\" not written" << std::endl;
    abort();
  }
  virtual double operator()(int c, const Vertex3d &v) const {
    std::cerr << "\"DomainRightHandSide::operator()\" not written" << std::endl;
    abort();
  }

  virtual void operator()(VectorIterator b, const TestFunction &N,
                          const Vertex2d &v) const {
    for (int c = 0; c < GetNcomp(); c++) {
      b[c] += N.m() * (*this)(c, v);
    }
  }
  virtual void operator()(VectorIterator b, const TestFunction &N,
                          const Vertex3d &v) const {
    for (int c = 0; c < GetNcomp(); c++) {
      b[c] += N.m() * (*this)(c, v);
    }
  }

  virtual void SetCellSize(double h) const {}
  virtual void point_cell(int material) const {}
};

typedef DomainRightHandSide DomainInitialCondition;

/*-------------------------------------------------------*/

} // namespace Gascoigne

#endif
