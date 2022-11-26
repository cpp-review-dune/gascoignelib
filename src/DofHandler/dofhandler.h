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

#include "../Common/paramfile.h"
#include "../Common/vertex.h"
#include "../Interface/gascoigne.h"
#include "../Mesh/boundaryindexhandler.h"
#include "../Mesh/hangingindexhandler.h"
#include "../Mesh/patchindexhandler.h"

/*-----------------------------------------*/

namespace Gascoigne {
class DofHandlerBase
{
protected:
  // VANKA: muss weg
  IndexVector mat_Vanka, mat_Vanka_patch;
  std::vector<std::array<Vertex3d, 3>> basis_Vanka, basis_Vanka_patch;

  IndexVector nc, vo2n, mat, matpatch;
  PatchIndexHandler PatchHandler;
  BoundaryIndexHandler BoundaryHandler;
  HangingIndexHandler HangingHandler;

public:
  //////////////////////////////////////////////////
  // Constructor, Init
  DofHandlerBase(){};
  virtual ~DofHandlerBase() {}

  virtual void BasicInit(const ParamFile* /*pf*/)
  {
    std::cerr << "\"DofHandler:BasicInit\" not written!" << std::endl;
    abort();
  }

  //////////////////////////////////////////////////
  // Access
  virtual IndexType dimension() const = 0;
  virtual IndexType nnodes() const = 0;
  virtual IndexType nhanging() const { return 0; }
  virtual IndexType ncells() const = 0;
  virtual IndexType nelements(IndexType degree) const = 0;

  virtual IndexType nodes_per_cell(IndexType i) const = 0;
  virtual IndexType vertex_of_cell(IndexType i, IndexType ii) const = 0;
  virtual const Vertex2d& vertex2d(IndexType /*i*/) const
  {
    std::cerr << "\"MeshInterface::vertex2d\" not written!" << std::endl;
    abort();
  }
  virtual const Vertex3d& vertex3d(IndexType /*i*/) const
  {
    std::cerr << "\"MeshInterface::vertex3d\" not written!" << std::endl;
    abort();
  }

  virtual const IndexVector& GetCellVector() const { return nc; }
  virtual const IndexVector& GetMaterialVector() const { return mat; }
  virtual const IndexVector& GetMaterialPatchVector() const { return matpatch; }

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
  virtual const IndexVector* Vertexo2n() const { return &vo2n; }

  virtual IndexVector& GetCellVector() { return nc; }
  virtual IndexVector& GetMaterialVector() { return mat; }
  virtual IndexVector& GetMaterialPatchVector() { return matpatch; }

  virtual PatchIndexHandler& GetPatchIndexHandler() { return PatchHandler; }
  virtual BoundaryIndexHandler& GetBoundaryIndexHandler()
  {
    return BoundaryHandler;
  }
  virtual HangingIndexHandler& GetHangingIndexHandler()
  {
    return HangingHandler;
  }
  virtual IndexVector* Vertexo2n() { return &vo2n; }
  virtual IndexType VtkType(IndexType i) const = 0;

  // wird von DofHandler2d/3d ueberschrieben
  virtual IndexVector IndicesOfCell(IndexType /*iq*/) const
  {
    std::cerr << "\"DofHandler:IndicesOfCell\" not written!" << std::endl;
    abort();
  }

  virtual bool CellIsCurved(IndexType /*iq*/) const { return 0; }
  virtual IndexType nedges() const { return 0; }

  virtual IndexVector GetElement(IndexType /*degree*/, IndexType /*iq*/) const
  {
    std::cerr << "DofHandler::GetElement" << std::endl;
    abort();
  }
  virtual IndexType nodes_per_element(IndexType /*degree*/) const
  {
    std::cerr << "DofHandler::nodes_per_element" << std::endl;
    abort();
  }

  // Access: Patch-Structure
  virtual IndexType nodes_per_patch() const
  {
    return PatchHandler.nodes_per_patch();
  }
  virtual IndexType nodes_per_q4patch() const
  {
    return PatchHandler.nodes_per_q4patch();
  }
  virtual bool HasPatch() const { return PatchHandler.HasPatch(); }
  virtual bool HasQ4Patch() const { return PatchHandler.HasQ4Patch(); }
  virtual IndexType npatches() const { return PatchHandler.npatches(); }
  virtual IndexType nq4patches() const { return PatchHandler.nq4patches(); }
  virtual const IndexVector* IndicesOfPatch(IndexType i) const
  {
    return &(PatchHandler.IndicesOfPatch(i));
  }
  virtual const IndexVector* IndicesOfQ4Patch(IndexType i) const
  {
    return &(PatchHandler.IndicesOfQ4Patch(i));
  }
  virtual IndexVector Q2IndicesOfQ4Patch(IndexType i) const
  {
    return PatchHandler.Q2IndicesOfQ4Patch(i);
  }
  virtual IndexVector CoarseIndices(IndexType iq) const
  {
    return PatchHandler.CoarseIndices(iq);
  }
  virtual IndexVector CoarseIndicesQ4(IndexType iq) const
  {
    return PatchHandler.CoarseIndicesQ4(iq);
  }

  // Access: material
  virtual IndexType material(IndexType i) const
  {
    assert(i < mat.size());
    return mat[i];
  }

  virtual IndexType material_patch(IndexType i) const
  {
    assert(i < matpatch.size());
    return matpatch[i];
  }

  IndexType material(IndexType degree, IndexType i) const
  {
    if (degree == 1)
      return material(i);
    else if (degree == 2)
      return material_patch(i);
    abort();
  }

  virtual IndexType& material(IndexType i)
  {
    assert(i < mat.size());
    return mat[i];
  }

  virtual IndexType& material_patch(IndexType i)
  {
    assert(i < matpatch.size());
    return matpatch[i];
  }

  // Access: boundary
  virtual const IndexVector* VertexOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Verteces(color));
  }
  virtual const IndexVector* CellOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Cells(color));
  }
  virtual const IndexVector* LocalOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Localind(color));
  }
  virtual const IndexVector* PatchOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.Patches(color));
  }
  virtual const IndexVector* LocalPatchOnBoundary(IndexType color) const
  {
    return &(BoundaryHandler.LocalPatchind(color));
  }

  virtual const IndexVector* ElementOnBoundary(IndexType degree,
                                               IndexType color) const
  {
    if (degree == 1)
      return CellOnBoundary(color);
    else if (degree == 2)
      return PatchOnBoundary(color);
    else
      abort();
  }
  virtual const IndexVector* ElementLocalOnBoundary(IndexType degree,
                                                    IndexType color) const
  {
    if (degree == 1)
      return LocalOnBoundary(color);
    else if (degree == 2)
      return LocalPatchOnBoundary(color);
    else
      abort();
  }

  virtual IndexSet GetColors() const { return BoundaryHandler.GetColors(); }

  virtual IndexType& material_Vanka(IndexType i)
  {
    assert(i < mat_Vanka.size());
    return mat_Vanka[i];
  }
  virtual IndexType& material_Vanka_patch(IndexType i)
  {
    assert(i < mat_Vanka_patch.size());
    return mat_Vanka_patch[i];
  }
  virtual std::array<Vertex3d, 3>& Vanka_basis(IndexType i)
  {
    assert(i < basis_Vanka.size());
    return basis_Vanka[i];
  }
  virtual std::array<Vertex3d, 3>& Vanka_basis_patch(IndexType i)
  {
    assert(i < basis_Vanka_patch.size());
    return basis_Vanka[i];
  }

  // VANKA: muss weg
  virtual const IndexVector& GetMaterialVankaVector() const
  {
    return mat_Vanka;
  }
  virtual const IndexVector& GetMaterialVankaPatchVector() const
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
  virtual IndexVector& GetMaterialVankaVector() { return mat_Vanka; }
  virtual IndexVector& GetMaterialVankaPatchVector() { return mat_Vanka_patch; }
  virtual std::vector<std::array<Vertex3d, 3>>& GetbasisVankaVector()
  {
    return basis_Vanka;
  }
  virtual std::vector<std::array<Vertex3d, 3>>& GetbasisVankaPatchVector()
  {
    return basis_Vanka_patch;
  }

  virtual IndexType material_Vanka(IndexType i) const
  {
    assert(i < mat_Vanka.size());
    return mat_Vanka[i];
  }
  virtual IndexType material_Vanka_patch(IndexType i) const
  {
    assert(i < mat_Vanka_patch.size());
    return mat_Vanka_patch[i];
  }
  virtual const std::array<Vertex3d, 3> Vanka_basis(IndexType i) const
  {
    assert(i < basis_Vanka.size());
    return basis_Vanka[i];
  }
  const std::array<Vertex3d, 3> Vanka_basis_patch(IndexType i) const
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

  IndexType dimension() const { return DIM; }
  IndexType nnodes() const { return nx.size(); }
  IndexType ncells() const { return nc.size() / ((DIM == 2) ? 4 : 8); }
  IndexType nelements(IndexType degree) const
  {
    if (degree == 1)
      return ncells();
    else if (degree == 2)
      return npatches();
    else
      abort();
  }

  IndexType nhanging() const { return HangingHandler.GetStructure()->size(); }

  IndexType nodes_per_cell(IndexType /*i*/) const { return (DIM == 2) ? 4 : 8; }
  IndexType nodes_per_element(IndexType degree) const
  {
    if (degree == 1)
      return nodes_per_cell(0);
    else if (degree == 2)
      return nodes_per_patch();
    else
      abort();
  }

  IndexType VtkType(IndexType /*i*/) const { return (DIM == 2) ? 9 : 12; }

  const Vertex<DIM>& vertex(IndexType i) const
  {
    assert(i < nx.size());
    return nx[i];
  }

  virtual const Vertex<2>& vertex2d(IndexType i) const;
  virtual const Vertex<3>& vertex3d(IndexType i) const;

  IndexType vertex_of_cell(IndexType i, IndexType ii) const
  {
    assert(nodes_per_cell(i) * i + ii < nc.size());
    return nc[nodes_per_cell(i) * i + ii];
  }

  virtual IndexType CornerIndices(IndexType degree,
                                  IndexType iq,
                                  IndexType ii) const
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

  IndexVector IndicesOfCell(IndexType iq) const;

  IndexVector GetElement(IndexType degree, IndexType iq) const
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
