/*----------------------------   dofhandler.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dofhandler_H
#define __dofhandler_H
/*----------------------------   dofhandler.h     ---------------------------*/

/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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

#include "boundaryindexhandler.h"
#include "gascoigne.h"
#include "hangingindexhandler.h"
#include "paramfile.h"
#include "patchindexhandler.h"
#include "vertex.h"

/*-----------------------------------------*/

namespace Gascoigne {
class DofHandlerBase
{
protected:
  // VANKA: muss weg
  IntVector mat_Vanka, mat_Vanka_patch;
  std::vector<std::array<Vertex3d, 3>> basis_Vanka, basis_Vanka_patch;

  IntVector nc, vo2n, mat, matpatch;
  PatchIndexHandler PatchHandler;
  BoundaryIndexHandler BoundaryHandler;
  HangingIndexHandler HangingHandler;

public:
  //////////////////////////////////////////////////
  // Constructor, Init
  DofHandlerBase(){};
  virtual ~DofHandlerBase() {}

  virtual void BasicInit(const ParamFile* pf)
  {
    std::cerr << "\"DofHandler:BasicInit\" not written!" << std::endl;
    abort();
  }

  //////////////////////////////////////////////////
  // Access
  virtual int dimension() const = 0;
  virtual int nnodes() const = 0;
  virtual int nhanging() const { return 0; }
  virtual int ncells() const = 0;
  virtual int nelements(int degree) const = 0;

  virtual int nodes_per_cell(int i) const = 0;
  virtual int vertex_of_cell(int i, int ii) const = 0;
  virtual const Vertex2d& vertex2d(int i) const
  {
    std::cerr << "\"MeshInterface::vertex2d\" not written!" << std::endl;
    abort();
  }
  virtual const Vertex3d& vertex3d(int i) const
  {
    std::cerr << "\"MeshInterface::vertex3d\" not written!" << std::endl;
    abort();
  }

  virtual const IntVector& GetCellVector() const { return nc; }
  virtual const IntVector& GetMaterialVector() const { return mat; }
  virtual const IntVector& GetMaterialPatchVector() const { return matpatch; }

  virtual const PatchIndexHandler& GetPatchIndexHandler() const
  {
    return PatchHandler;
  }
  virtual const BoundaryIndexHandler& GetBoundaryIndexHandler() const
  {
    return BoundaryHandler;
  }
  virtual const HangingIndexHandler& GetHangingIndexHandler() const
  {
    return HangingHandler;
  }
  virtual const IntVector* Vertexo2n() const { return &vo2n; }

  virtual IntVector& GetCellVector() { return nc; }
  virtual IntVector& GetMaterialVector() { return mat; }
  virtual IntVector& GetMaterialPatchVector() { return matpatch; }

  virtual PatchIndexHandler& GetPatchIndexHandler() { return PatchHandler; }
  virtual BoundaryIndexHandler& GetBoundaryIndexHandler()
  {
    return BoundaryHandler;
  }
  virtual HangingIndexHandler& GetHangingIndexHandler()
  {
    return HangingHandler;
  }
  virtual IntVector* Vertexo2n() { return &vo2n; }
  virtual int VtkType(int i) const = 0;

  // wird von DofHandler2d/3d ueberschrieben
  virtual IntVector IndicesOfCell(int iq) const
  {
    std::cerr << "\"DofHandler:IndicesOfCell\" not written!" << std::endl;
    abort();
  }

  virtual bool CellIsCurved(int iq) const { return 0; }
  virtual int nedges() const { return 0; }

  virtual IntVector GetElement(int degree, int iq) const
  {
    std::cerr << "DofHandler::GetElement" << std::endl;
    abort();
  }
  virtual int nodes_per_element(int degree) const
  {
    std::cerr << "DofHandler::nodes_per_element" << std::endl;
    abort();
  }

  // Access: Patch-Structure
  virtual int nodes_per_patch() const { return PatchHandler.nodes_per_patch(); }
  virtual int nodes_per_q4patch() const
  {
    return PatchHandler.nodes_per_q4patch();
  }
  virtual bool HasPatch() const { return PatchHandler.HasPatch(); }
  virtual bool HasQ4Patch() const { return PatchHandler.HasQ4Patch(); }
  virtual int npatches() const { return PatchHandler.npatches(); }
  virtual int nq4patches() const { return PatchHandler.nq4patches(); }
  virtual const IntVector* IndicesOfPatch(int i) const
  {
    return &(PatchHandler.IndicesOfPatch(i));
  }
  virtual const IntVector* IndicesOfQ4Patch(int i) const
  {
    return &(PatchHandler.IndicesOfQ4Patch(i));
  }
  virtual IntVector Q2IndicesOfQ4Patch(int i) const
  {
    return PatchHandler.Q2IndicesOfQ4Patch(i);
  }
  virtual IntVector CoarseIndices(int iq) const
  {
    return PatchHandler.CoarseIndices(iq);
  }
  virtual IntVector CoarseIndicesQ4(int iq) const
  {
    return PatchHandler.CoarseIndicesQ4(iq);
  }

  // Access: material
  virtual const int material(int i) const
  {
    assert(i < mat.size());
    return mat[i];
  }
  virtual const int material_patch(int i) const
  {
    assert(i < matpatch.size());
    return matpatch[i];
  }
  const int material(int degree, int i) const
  {
    if (degree == 1)
      return material(i);
    else if (degree == 2)
      return material_patch(i);
    abort();
  }
  virtual int& material(int i)
  {
    assert(i < mat.size());
    return mat[i];
  }
  virtual int& material_patch(int i)
  {
    assert(i < matpatch.size());
    return matpatch[i];
  }

  // Access: boundary
  virtual const IntVector* VertexOnBoundary(int color) const
  {
    return &(BoundaryHandler.Verteces(color));
  }
  virtual const IntVector* CellOnBoundary(int color) const
  {
    return &(BoundaryHandler.Cells(color));
  }
  virtual const IntVector* LocalOnBoundary(int color) const
  {
    return &(BoundaryHandler.Localind(color));
  }
  virtual const IntVector* PatchOnBoundary(int color) const
  {
    return &(BoundaryHandler.Patches(color));
  }
  virtual const IntVector* LocalPatchOnBoundary(int color) const
  {
    return &(BoundaryHandler.LocalPatchind(color));
  }

  virtual const IntVector* ElementOnBoundary(int degree, int color) const
  {
    if (degree == 1)
      return CellOnBoundary(color);
    else if (degree == 2)
      return PatchOnBoundary(color);
    else
      abort();
  }
  virtual const IntVector* ElementLocalOnBoundary(int degree, int color) const
  {
    if (degree == 1)
      return LocalOnBoundary(color);
    else if (degree == 2)
      return LocalPatchOnBoundary(color);
    else
      abort();
  }

  virtual std::set<int> GetColors() const
  {
    return BoundaryHandler.GetColors();
  }

  virtual int& material_Vanka(int i)
  {
    assert(i < mat_Vanka.size());
    return mat_Vanka[i];
  }
  virtual int& material_Vanka_patch(int i)
  {
    assert(i < mat_Vanka_patch.size());
    return mat_Vanka_patch[i];
  }
  virtual std::array<Vertex3d, 3>& Vanka_basis(int i)
  {
    assert(i < basis_Vanka.size());
    return basis_Vanka[i];
  }
  virtual std::array<Vertex3d, 3>& Vanka_basis_patch(int i)
  {
    assert(i < basis_Vanka_patch.size());
    return basis_Vanka[i];
  }

  // VANKA: muss weg
  virtual const IntVector& GetMaterialVankaVector() const { return mat_Vanka; }
  virtual const IntVector& GetMaterialVankaPatchVector() const
  {
    return mat_Vanka_patch;
  }
  virtual const std::vector<std::array<Vertex3d, 3>>& GetbasisVankaVector()
    const
  {
    return basis_Vanka;
  }
  virtual const std::vector<std::array<Vertex3d, 3>>& GetbasisVankaPatchVector()
    const
  {
    return basis_Vanka_patch;
  }
  virtual IntVector& GetMaterialVankaVector() { return mat_Vanka; }
  virtual IntVector& GetMaterialVankaPatchVector() { return mat_Vanka_patch; }
  virtual std::vector<std::array<Vertex3d, 3>>& GetbasisVankaVector()
  {
    return basis_Vanka;
  }
  virtual std::vector<std::array<Vertex3d, 3>>& GetbasisVankaPatchVector()
  {
    return basis_Vanka_patch;
  }

  virtual const int material_Vanka(int i) const
  {
    assert(i < mat_Vanka.size());
    return mat_Vanka[i];
  }
  virtual const int material_Vanka_patch(int i) const
  {
    assert(i < mat_Vanka_patch.size());
    return mat_Vanka_patch[i];
  }
  virtual const std::array<Vertex3d, 3> Vanka_basis(int i) const
  {
    assert(i < basis_Vanka.size());
    return basis_Vanka[i];
  }
  const std::array<Vertex3d, 3> Vanka_basis_patch(int i) const
  {
    assert(i < basis_Vanka_patch.size());
    return basis_Vanka_patch[i];
  }
};

template<int DIM>
class DofHandler : public DofHandlerBase
{
protected:
  // basic
  std::vector<Vertex<DIM>> nx;

public:
  DofHandler(){};
  ~DofHandler(){};

  std::string GetName() const;
  std::vector<Vertex<DIM>>& GetVertexVector() { return nx; }
  const std::vector<Vertex<DIM>>& GetVertexVector() const { return nx; }

  int dimension() const { return DIM; }
  int nnodes() const { return nx.size(); }
  int ncells() const { return nc.size() / ((DIM == 2) ? 4 : 8); }
  int nelements(int degree) const
  {
    if (degree == 1)
      return ncells();
    else if (degree == 2)
      return npatches();
    else
      abort();
  }

  int nhanging() const { return HangingHandler.GetStructure()->size(); }

  int nodes_per_cell(int i) const { return (DIM == 2) ? 4 : 8; }
  int nodes_per_element(int degree) const
  {
    if (degree == 1)
      return nodes_per_cell(0);
    else if (degree == 2)
      return nodes_per_patch();
    else
      abort();
  }

  int VtkType(int i) const { return (DIM == 2) ? 9 : 12; }

  const Vertex<DIM>& vertex(int i) const
  {
    assert(i < nx.size());
    return nx[i];
  }

  virtual const Vertex<2>& vertex2d(int i) const;
  virtual const Vertex<3>& vertex3d(int i) const;

  int vertex_of_cell(int i, int ii) const
  {
    assert(nodes_per_cell(i) * i + ii < nc.size());
    return nc[nodes_per_cell(i) * i + ii];
  }

  virtual int CornerIndices(int degree, int iq, int ii) const
  {
    if (degree == 1) {
      assert(nodes_per_cell(iq) * iq + ii < nc.size());
      return nc[nodes_per_cell(iq) * iq + ii];
    } else if (degree == 2) {
      return PatchHandler.CoarseIndices(iq)[ii];
    }
    assert(0);
    return 0;
  }

  IntVector IndicesOfCell(int iq) const;

  IntVector GetElement(int degree, int iq) const
  {
    if (degree == 1)
      return IndicesOfCell(iq);
    else if (degree == 2)
      return *IndicesOfPatch(iq);
    abort();
  }
};

} // namespace Gascoigne

/*----------------------------   dofhandler.h     ---------------------------*/
/* end of #ifndef __dofhandler_H */
#endif
/*----------------------------   dofhandler.h     ---------------------------*/
