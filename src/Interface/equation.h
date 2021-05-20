/**
 *
 * Copyright (C) 2004, 2005,2018 by the Gascoigne 3D authors
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

#ifndef __Equation_h
#define __Equation_h

#include "application.h"
#include "entrymatrix.h"
#include "compvector.h"
#include "vertex.h"

/*-------------------------------------------------------------------------*/

namespace Gascoigne {

/**
  \brief Interface class for defining equations

  Defines an equation object in Gascoigne. Equations may not have persistent
  data, only references or pointers to a data object are ok.

  The most important members are

  - Form(U,N)
    which defines the variational formulation, U being the trial-functions, N
  the tests-functions. Form fills a vector b of size ncomp (number of solutions
  components). Remember to allways add values to b[...]

  - Matrix(U,M,N)
    defines the Jacobian of the variational formulation at U. Fills a matrix A
  of size ncomp x ncomps. A(j,i) corresponds to test-function number i and
  derivative in direction U[j]

  - point(U,v)
    for precomputing values to be used in Form. Depends on solution state U and
  coordinate v, but is reused for all test-functions

  - pointmatrix(U,v)
    same but used in Matrix. If pointmatrix() is not given, point() is used
  instead

  - createNew()
    used to clone the entire class. Required for parallelization where each
  trhead has it own instance of the equation object.

 */

class Equation : public virtual Application {
private:
protected:
public:
  //
  // Constructors
  //
  Equation() : Application() {}
  virtual ~Equation() {}

  /**
     clones an Equation. Usually it is simply to return a new instance of the
     same object passing the required variables to the constructor. It takes the
     role of a copy constructor and the cloning of classes is required for
     multithreading.
  */
  virtual Equation *createNew() const {
    std::cerr << "\"Equation::createNew\" not written!" << std::endl;
    abort();
  }

  virtual void OperatorStrong(DoubleVector &b, const FemFunction &U) const {
    std::cerr << "\"Equation::OperatorStrong\" not written!" << std::endl;
    abort();
  }
  virtual void Pattern(TimePattern &TP) const {
    std::cerr << "\"Equation::SetTimePattern\" not written!" << std::endl;
    abort();
  }

  virtual void point_cell(int material) const {}

  /// point is called once the quadrature point is defined and the finite
  /// element function U is initialized in this point.
  virtual void point(double h, const FemFunction &U, const Vertex2d &v) const {}
  virtual void point(double h, const FemFunction &U, const Vertex3d &v) const {}

  /// point_M is after calling point(...) and before calling Matrix(...) here,
  /// the update M is initialized in the quadrature point and more values can
  /// be precomputed
  virtual void point_M(const FemFunction &U, const TestFunction &M) const {}

  virtual void pointmatrix(double h, const FemFunction &U,
                           const Vertex2d &v) const {
    point(h, U, v);
  }
  virtual void pointmatrix(double h, const FemFunction &U,
                           const Vertex3d &v) const {
    point(h, U, v);
  }

  // Elements
  virtual int GetNcomp() const = 0;
  virtual void Form(VectorIterator b, const FemFunction &U,
                    const TestFunction &N) const = 0;
  virtual void Matrix(EntryMatrix &A, const FemFunction &U,
                      const TestFunction &M, const TestFunction &N) const = 0;

  // Boundary
  virtual void Form(VectorIterator b, const FemFunction &U,
                    const TestFunction &N, int col) const {}
  virtual void Matrix(EntryMatrix &E, const FemFunction &U,
                      const TestFunction &M, const TestFunction &N,
                      int col) const {}

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

  virtual void MatrixBlock(EntryMatrix &A, const FemFunction &U,
                           const FemFunction &N) const {
    for (int j = 0; j < N.size(); j++) {
      point_M(U, N[j]);
      for (int i = 0; i < N.size(); i++) {
        A.SetDofIndex(i, j);
        Matrix(A, U, N[j], N[i]);
      }
    }
  }
};
} // namespace Gascoigne

#endif
