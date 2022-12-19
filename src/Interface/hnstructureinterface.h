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

#ifndef __HNStructureInterface_h
#define __HNStructureInterface_h

#include "../DofHandler/dofhandlerbase.h"
#include "../LinAlg/sparsestructure.h"

#include "gascoigne.h"
#include "matrixinterface.h"

/*-----------------------------------------*/

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments HNStructureInterface

////
////
/////////////////////////////////////////////

class HNStructureInterface
{
private:
protected:
public:
  //
  ////  Con(De)structor
  //
  HNStructureInterface() {}
  virtual ~HNStructureInterface() {}

  virtual void ReInit(const DofHandlerBase* m) = 0;
  virtual void MatrixDiag(ShortIndexType ncomp, MatrixInterface& A) const = 0;
  virtual void SparseStructureDiag(SparseStructure* S) const = 0;

  virtual void CondenseHanging(IndexVector& indices) const = 0;
  virtual void CondenseHanging(EntryMatrix& E, IndexVector& indices) const = 0;
  virtual void CondenseHangingPatch(EntryMatrix& E,
                                    IndexVector& indices) const = 0;

  virtual void Average(GlobalVector& u) const = 0;
  virtual void Distribute(GlobalVector& u) const = 0;
  virtual void Zero(GlobalVector& u) const = 0;
  virtual bool ZeroCheck(const GlobalVector& u) const = 0;

  virtual IndexType nhnodes() const = 0;
};
} // namespace Gascoigne

#endif
