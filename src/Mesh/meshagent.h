/**
 *
 * Copyright (C) 2004, 2005, 2006, 2009, 2011 by the Gascoigne 3D authors
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

#ifndef __MeshAgent_h
#define __MeshAgent_h

#include "boundaryfunction.h"
#include "gascoignemultigridmesh.h"
#include "hierarchicalmesh.h"
#include "stdperiodicmapping.h"

/*-----------------------------------------*/

namespace Gascoigne {
typedef std::map<int, int> IntMap;

class MeshAgent
{
private:
  std::map<int, BoundaryFunction<2>*> _curved2d;
  std::map<int, BoundaryFunction<3>*> _curved3d;

  // Fuer die Zuordnung GM Nr auf altem Gitter zu GM Nr. auf neuem Gitter
  IntVector _cl2g, _celll2g;
  IntVector _fathers; // GM Nr zu HM nr.
  IntMap _cg2l, _cellg2l;
  std::map<int, std::set<int>> _co2n;
  bool _goc2nc;
  nvector<IntVector> _q4patch, _q4toq2;

protected:
  HierarchicalMesh* HMP;
  GascoigneMultiGridMesh* GMG;
  IntVector _periodicCols;
  std::map<int, std::map<int, PeriodicMapping*>> _periodicMaps;

public:
  MeshAgent();
  virtual ~MeshAgent();

  virtual void ReInit();

  virtual int GetDimension() const { return HMP->dimension(); }

  virtual GascoigneMultiGridMesh* NewMultiGridMesh()
  {
    return new GascoigneMultiGridMesh;
  }

  virtual void BuildQ4PatchList(const IntVector& patchl2g);

  virtual void AssemblePeriodicBoundaries();

  virtual GascoigneMesh* GMesh(int l) { return GMG->GetGascoigneMesh(l); }

  virtual void AddShape(int col, BoundaryFunction<2>* f) { _curved2d[col] = f; }
  virtual void AddShape(int col, BoundaryFunction<3>* f) { _curved3d[col] = f; }

  virtual std::map<int, BoundaryFunction<2>*>& GetShapes2d()
  {
    return _curved2d;
  }
  virtual std::map<int, BoundaryFunction<3>*>& GetShapes3d()
  {
    return _curved3d;
  }
  virtual const std::map<int, BoundaryFunction<2>*>& GetShapes2d() const
  {
    return _curved2d;
  }
  virtual const std::map<int, BoundaryFunction<3>*>& GetShapes3d() const
  {
    return _curved3d;
  }

  virtual void AddPeriodicMapping(int col, int col2, PeriodicMapping* map)
  {
    _periodicMaps[col][col2] = map;
  }

  virtual void BasicInit(const ParamFile& pf);
  virtual void BasicInit(const ParamFile& pf, int pdepth);
  virtual void BasicInit(const std::string& gridname,
                         int dim,
                         int patchdepth,
                         int epatcher,
                         bool goc2nc = false);

  virtual const GascoigneMultiGridMesh& GetMultiGrid() const { return *GMG; }
  virtual GascoigneMultiGridMesh& GetMultiGrid() { return *GMG; }

  virtual HierarchicalMesh* GetHierarchicalMesh() { return HMP; }
  virtual const HierarchicalMesh* GetHierarchicalMesh() const { return HMP; }

  virtual int nnodes() const { return GMG->GetGascoigneMesh(0)->nnodes(); }
  virtual int ncells() const { return GMG->GetGascoigneMesh(0)->ncells(); }
  virtual int nlevels() const { return GMG->nlevels(); }

  virtual const GascoigneMesh* GetMesh() const
  {
    return GMG->GetGascoigneMesh(0);
  }
  virtual const GascoigneMesh* GetMesh(int l) const
  {
    return GMG->GetGascoigneMesh(l);
  }

  virtual void read_gup(const std::string& fname);
  virtual void read_gip(const std::string& fname);
  virtual void write_gup(const std::string& fname) const;
  virtual void write_gip(const std::string& fname) const;
  virtual void write_inp(const std::string& fname) const;
  virtual void global_refine(int n);
  virtual void global_patch_coarsen(int n);
  virtual void random_patch_coarsen(double p, int n);
  virtual void random_patch_refine(double p, int n);
  virtual void refine_nodes(IntVector& refnodes, IntVector& coarsenodes);
  virtual void refine_nodes(IntVector& refnodes);
  virtual void refine_cells(IntVector& ref);

  virtual const GascoigneMeshTransfer* GetTransfer(int l) const
  {
    return GMG->GetTransfer(l);
  }

  virtual const std::set<int> Cello2n(int i) const;
  virtual const int Cello2nFather(int i) const;
  virtual void ClearCl2g() { _cl2g.clear(); }
  virtual const bool Goc2nc() const { return _goc2nc; }

  virtual const IntVector& Celll2g() const { return _celll2g; }
  virtual const IntMap& Cellg2l() const { return _cellg2l; }
};
} // namespace Gascoigne

#endif
