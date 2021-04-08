/*----------------------------   laplace.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __laplace_H
#define __laplace_H
/*----------------------------   laplace.h     ---------------------------*/

/**
 *
 * Copyright (C) 2004,2018 by the Gascoigne 3D authors
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

#include "equation.h"
#include "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne {

class LaplaceData {
public:
  double visc;

  void BasicInit(const ParamFile &pf) {
    DataFormatHandler DFH;
    DFH.insert("visc", &visc, 1.);
    FileScanner FS(DFH, pf, "Equation");
  }
};

template <int DIM> class Laplace : public virtual Equation {

protected:
  LaplaceData data;

public:
  Laplace(const LaplaceData &PD) : data(PD) {}

  Laplace<DIM> *createNew() const { return new Laplace<DIM>(data); }

  std::string GetName() const { return "Laplace"; }
  int GetNcomp() const { return 1; }
  void point(double h, const FemFunction &U, const Vertex2d &v) const {}

  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i)
      b[0] += data.visc * U[0][i + 1] * N[i + 1];
  }

  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const {
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += data.visc * M[i + 1] * N[i + 1];
  }
};

#define Laplace3d Laplace<3>
#define Laplace2d Laplace<2>
} // namespace Gascoigne

/*----------------------------   laplace.h     ---------------------------*/
/* end of #ifndef __laplace_H */
#endif
/*----------------------------   laplace.h     ---------------------------*/
