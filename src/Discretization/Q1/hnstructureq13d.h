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

#ifndef __HNStructureQ13d_h
#define __HNStructureQ13d_h

#include "../../Common/entrymatrix.h"

#include "hnstructureq12d.h"

/*-----------------------------------------*/

namespace Gascoigne {
class HNStructureQ13d : public HNStructureQ12d
{
protected:
  typedef std::array<IndexType, 9> FaceVector;

  typedef std::map<IndexType, FaceVector>::iterator fiterator;
  typedef std::map<IndexType, FaceVector>::const_iterator const_fiterator;

  const std::map<IndexType, FaceVector>* faces;

  std::array<std::array<IndexType, 3>, 12> lnoe;
  std::array<std::array<IndexType, 5>, 6> lnop;

  void CondenseHanging2er(IndexVector& indices) const;
  void CondenseHanging4er(IndexVector& indices) const;

  void CondenseHanging2er(EntryMatrix& E, IndexVector& indices) const;
  void CondenseHanging4er(EntryMatrix& E, IndexVector& indices) const;

  std::array<IndexType, 4> GetHangingFace(IndexType i) const;
  std::array<IndexType, 2> GetHangingEdge(IndexType i) const;

public:
  ~HNStructureQ13d();
  HNStructureQ13d();

  IndexType nhnodes() const { return edges->size() + faces->size(); }

  void ReInit(const GascoigneMesh* m);
  IndexType hanging(IndexType i) const;

  void MatrixDiag(ShortIndexType ncomp, MatrixInterface& A) const;
  void SparseStructureDiag(SparseStructure* A) const;

  void Average(GlobalVector& u) const;
  void Distribute(GlobalVector& u) const;
  void Zero(GlobalVector& u) const;
  bool ZeroCheck(const GlobalVector& u) const;

  void Couplings(IndexVector& indices) const;

  void CondenseHanging(IndexVector& indices) const;
  void CondenseHanging(EntryMatrix&, IndexVector&) const;
  void CondenseHangingPatch(EntryMatrix& E, IndexVector& indices) const;
};
} // namespace Gascoigne

#endif
