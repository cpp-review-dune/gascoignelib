/**
 *
 * Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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

#ifndef __p4estmesh2d_h
#define __p4estmesh2d_h

#include <set>
#include <string>
#include <utility>

#include <p4est_bits.h>
#include <p4est_vtk.h>

#include "../Common/compvector.h"
#include "../Common/paramfile.h"
#include "../Common/triple.h"
#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"

#include "curvedshapes.h"
#include "edge.h"
#include "hanglist.h"
#include "hierarchicalmesh2d.h"

/*---------------------------------------------------*/

namespace Gascoigne {
class P4estMesh2d : public HierarchicalMesh2d
{
protected:
  /*  typedef  */

  typedef std::pair<int, int> pint;
  typedef triple<int, int, int> tint;
  typedef std::array<int, 2> EdgeVector;
  typedef std::array<int, 4> FaceVector;
  typedef std::vector<Edge> EdgeVec;
  typedef IntSet::iterator IntSetIt;
  typedef IntSet::const_iterator IntSetCIt;

  /*  Data  */

  p4est_t* p4est;
  p4est_connectivity_t* conn;

  void update_edges(IntVector&) { NOT_IMPLEMENTED }
  int FindPatchDepth() const { NOT_IMPLEMENTED };
  void FillVertexLevels(IntVector& dst) const { NOT_IMPLEMENTED };
  void RefineCoarseNodes(IntSet& dst,
                         const IntVector& refnodes,
                         const IntVector& vertexlevel) const {
    NOT_IMPLEMENTED
  };
  void VertexToCells(IntVector& dst,
                     const IntSet& src,
                     const IntVector& vertexlevel) const { NOT_IMPLEMENTED };
  void VertexToCellsCoarsening(IntVector& dst,
                               const IntSet& src,
                               const IntVector& vertexlevel) const {
    NOT_IMPLEMENTED
  };

public:
  ~P4estMesh2d()
  {
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
  };

  int withfaces;

  P4estMesh2d(){};
  P4estMesh2d(const P4estMesh2d&){ NOT_IMPLEMENTED };
  P4estMesh2d& operator=(const P4estMesh2d&){ NOT_IMPLEMENTED };

  /*  Zugriff  */

  int nnodes() const { NOT_IMPLEMENTED };
  int nlevels() const { NOT_IMPLEMENTED }
  int nedges() const { NOT_IMPLEMENTED }

  const IntVector* Vertexo2n() const { NOT_IMPLEMENTED }
  const IntVector* Edgeo2n() const { NOT_IMPLEMENTED }
  const IntVector* Cello2n() const { NOT_IMPLEMENTED }

  int Vertexo2n(int i) const { NOT_IMPLEMENTED }
  int Edgeo2n(int i) const { NOT_IMPLEMENTED }
  int Cello2n(int i) const { NOT_IMPLEMENTED }

  const Edge& edge(int i) const { NOT_IMPLEMENTED }
  const EdgeVec& edge() const { NOT_IMPLEMENTED }

  int child(int i, int ii) const { NOT_IMPLEMENTED };
  int nchilds(int i) const { NOT_IMPLEMENTED };

  int level(int i) const { NOT_IMPLEMENTED };
  bool sleep(int i) const { NOT_IMPLEMENTED };

  int Vater(const int i) const { NOT_IMPLEMENTED } IntVector
    Nachkommen(const int i) const { NOT_IMPLEMENTED } IntVector
    Geschwister(const int i) const { NOT_IMPLEMENTED } IntVector
    Kinder(const int i) const
  {
    NOT_IMPLEMENTED
  }

  int nodes_per_cell(int i) const { NOT_IMPLEMENTED };
  int vertex_of_cell(int i, int ii) const { NOT_IMPLEMENTED };

  void SetParameters(std::string gridname,
                     int patchdepth,
                     int epatcher){ NOT_IMPLEMENTED };
  void SetParameters(int patchdepth){ NOT_IMPLEMENTED };
  void ReadFile(const std::string& gridname){ NOT_IMPLEMENTED };
  void BasicInit(const ParamFile& pf, int pdepth = 0){ NOT_IMPLEMENTED };
  void global_refine(int k){ NOT_IMPLEMENTED };
  void global_patch_coarsen(int k){ NOT_IMPLEMENTED };
  void random_refine(double, int k = 1){ NOT_IMPLEMENTED };
  void random_patch_refine(double, int k = 1){ NOT_IMPLEMENTED };
  void random_patch_coarsen(double, int k = 0){ NOT_IMPLEMENTED };
  void random_double_patch_refine(double, int k = 1){ NOT_IMPLEMENTED };
  void clear_transfer_lists(){ NOT_IMPLEMENTED };
  void write_gip(const std::string&) const { NOT_IMPLEMENTED };
  void write_gup(const std::string&) const { NOT_IMPLEMENTED };
  void write_inp(const std::string&) const { NOT_IMPLEMENTED };

  int dimension() const { return 2; }
  int ncells() const { return p4est->global_num_quadrants; };
  int patchdepth() const { NOT_IMPLEMENTED }
  int nactivedescendants(int i) const { NOT_IMPLEMENTED };
  IntVector GetVertices(int c) const { NOT_IMPLEMENTED };

  bool CellIsCurved(int iq) const { NOT_IMPLEMENTED }

  std::set<int> GetColors() const { NOT_IMPLEMENTED };

  void read_inp(const std::string&); // IMPLEMENTED
  void read_gup(const std::string&){ NOT_IMPLEMENTED };
  void read_gip(const std::string&){ NOT_IMPLEMENTED };
  void refine(const IntVector&, const IntVector&); // IMPLEMENTED
  void patch_refine(IntVector&, IntVector&){ NOT_IMPLEMENTED };
  void vertex_patch_refine(IntVector& ref,
                           IntVector& coarse){ NOT_IMPLEMENTED };
  void vertex_patch_refine(IntVector&){ NOT_IMPLEMENTED };
  void GetAwakePatchs(std::set<int>&) const { NOT_IMPLEMENTED };
  void GetAwakeCells(std::set<int>&) const { TO_DO };
  void ConstructQ2PatchMesh(IntVector& pm) const { NOT_IMPLEMENTED };
  IntVector ConstructQ4Patch(int c) const { NOT_IMPLEMENTED };
  std::set<int> CellNeighbours(int i) const;

  int GetBoundaryCellOfCurved(int iq) const { NOT_IMPLEMENTED }

  void AddShape(int col, BoundaryFunction<2>* f) { NOT_IMPLEMENTED }
  void AddShape(int col, BoundaryFunction<3>* f) { NOT_IMPLEMENTED }

  void ShowOutput(int i) const { NOT_IMPLEMENTED }
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif //__p4estmesh2d_h
