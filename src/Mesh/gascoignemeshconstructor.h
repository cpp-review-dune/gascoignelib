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

#ifndef __GascoigneMeshConstructor_h
#define __GascoigneMeshConstructor_h

#include <map>

#include "../Interface/gascoigne.h"

#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "gascoignemultigridmesh.h"
#include "levelmesh2d.h"
#include "levelmesh3d.h"

/*-----------------------------------------*/

namespace Gascoigne {

class GascoigneMeshConstructor
{
private:
  IndexVector _cl2g, _pl2g;
  IndexMap _cg2l;

protected:
  const HierarchicalMesh* HM;
  GascoigneMultiGridMesh* GMG;

  bool finestlevel;

  virtual void PatchToCell2d(PatchIndexHandler& PIH,
                             const LevelMesh2d* LM) const;
  virtual void PatchToCell3d(PatchIndexHandler& PIH,
                             const LevelMesh3d* LM) const;

  virtual void Construct2d(GascoigneMesh* NM, const LevelMesh2d* LM) const;
  virtual void Construct3d(GascoigneMesh* NM, const LevelMesh3d* LM) const;

  virtual LevelMesh2d* LevelUpdate2d(GascoigneMesh* GM,
                                     const IndexSet& newquads,
                                     const IndexSet& oldquads) const;
  virtual LevelMesh3d* LevelUpdate3d(GascoigneMesh* GM,
                                     const IndexSet& newquads,
                                     const IndexSet& oldquads) const;
  virtual void Loop2d();
  virtual void Loop3d();

public:
  GascoigneMeshConstructor(const HierarchicalMesh* mm,
                           GascoigneMultiGridMesh* gmg);
  virtual ~GascoigneMeshConstructor() {}

  virtual void BasicInit();
  const IndexVector& Patchl2g() const { return _pl2g; }
  const IndexVector& Celll2g() const { return _cl2g; }
  const IndexMap& Cellg2l() const { return _cg2l; }
};
} // namespace Gascoigne

#endif
