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

#ifndef __MgInterpolatorMatrix_h
#define __MgInterpolatorMatrix_h

#include "columnstencil.h"
#include "gascoigne.h"
#include "meshtransferinterface.h"
#include "mginterpolatorinterface.h"

/*-----------------------------------------*/

namespace Gascoigne {
class MgInterpolatorMatrix : public virtual MgInterpolatorInterface
{
private:
  // We store two matrices, one for restrict, one for prolongate
  // to allow for parallelization of the matrix-vector product
  ColumnStencil STfine, STcoarse;
  DoubleVector valfine, valcoarse;

public:
  MgInterpolatorMatrix()
    : MgInterpolatorInterface()
  {
  }
  void BasicInit(const MeshTransferInterface* MT);

  void restrict_zero(GlobalVector& uL, const GlobalVector& ul) const;
  void prolongate_add(GlobalVector& ul, const GlobalVector& uL) const;
  void SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const;
};
} // namespace Gascoigne

#endif
