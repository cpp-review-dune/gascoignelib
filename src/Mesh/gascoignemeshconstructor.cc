/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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

#include "gascoignemeshconstructor.h"
#include "gascoignemeshtransferconstructor.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {

GascoigneMeshConstructor::GascoigneMeshConstructor(const HierarchicalMesh* hm,
                                                   GascoigneMultiGridMesh* gmg)
  : HM(hm)
  , GMG(gmg)
  , finestlevel(0)
{}

/*---------------------------------------------------*/

void
GascoigneMeshConstructor::BasicInit()
{
  if (HM->dimension() == 2)
    Loop2d();
  else
    Loop3d();

  IndexVector& v1 = *GMG->GetGascoigneMesh(0)->Vertexo2n();
  v1.reservesize(HM->Vertexo2n()->size());
  v1 = *HM->Vertexo2n();
}

/*---------------------------------------------------*/

void
GascoigneMeshConstructor::Loop2d()
{
  finestlevel = 1;

  IndexSet newquads, oldquads;
  HM->GetAwakeCells(newquads);

  GascoigneMesh* GM0 = GMG->GetGascoigneMesh(0);

  LevelMesh2d* FM = NULL;
  LevelMesh2d* CM = NULL;

  FM = LevelUpdate2d(GM0, newquads, oldquads);

  // Q4-PatchStruktur
  if (HM->patchdepth() >= 2) {
    PatchIndexHandler& PIH = GM0->GetPatchIndexHandler();
    FM->ConstructCellIndOfPatch(_pl2g);
    PIH.GetHasQ4Patch() = 1;
    HangingIndexHandler& H = GM0->GetHangingIndexHandler();
    FM->ConstructHangingStructureQuartic(*(H.GetQ4Structure()));
  }

  // Wir brauchen local2global und global2local Vector aus FM
  _cl2g = FM->Quadl2g();
  _cg2l = FM->Quadg2l();
  const HierarchicalMesh2d* HMP = dynamic_cast<const HierarchicalMesh2d*>(HM);
  assert(HMP);

  for (IndexType level = 1; level < GMG->nlevels(); level++) {
    finestlevel = 0;

    GascoigneMesh* GM = GMG->GetGascoigneMesh(level);
    GascoigneMeshTransfer* T = GMG->GetTransfer(level - 1);

    FM->construct_lists(newquads, oldquads);
    CM = LevelUpdate2d(GM, newquads, oldquads);

    GascoigneMeshTransferConstructor2d(HMP, T, FM, CM);

    delete FM;
    FM = CM;
    CM = NULL;
  }
  delete FM;
  FM = NULL;
}

/*---------------------------------------------------*/

void
GascoigneMeshConstructor::Loop3d()
{
  finestlevel = 1;

  IndexSet newquads, oldquads;
  HM->GetAwakeCells(newquads);

  GascoigneMesh* GM0 = GMG->GetGascoigneMesh(0);

  LevelMesh3d* FM = NULL;
  LevelMesh3d* CM = NULL;

  FM = LevelUpdate3d(GM0, newquads, oldquads);

  // Q4-PatchStruktur
  if (HM->patchdepth() >= 2) {
    PatchIndexHandler& PIH = GM0->GetPatchIndexHandler();
    FM->ConstructCellIndOfPatch(_pl2g);
    PIH.GetHasQ4Patch() = 1;
    HangingIndexHandler& H = GM0->GetHangingIndexHandler();
    FM->ConstructHangingStructureQuartic(*(H.GetQ4Structure()),
                                         *(H.GetQ4StructureFace()));
  }

  // Wir brauchen local2global und global2local Vector aus FM
  _cl2g = FM->Hexl2g();
  _cg2l = FM->Hexg2l();
  const HierarchicalMesh3d* HMP = dynamic_cast<const HierarchicalMesh3d*>(HM);
  assert(HMP);

  for (IndexType level = 1; level < GMG->nlevels(); level++) {
    finestlevel = 0;

    GascoigneMesh* GM = GMG->GetGascoigneMesh(level);
    GascoigneMeshTransfer* T = GMG->GetTransfer(level - 1);

    FM->construct_lists(newquads, oldquads);
    //      cout << "*** construct_lists " << level << "\t" << newquads.size()
    //      << " " << oldquads.size() << endl;

    CM = LevelUpdate3d(GM, newquads, oldquads);

    GascoigneMeshTransferConstructor3d(HMP, T, FM, CM);

    delete FM;
    FM = CM;
    CM = NULL;
  }
  delete FM;
  FM = NULL;
}

/*---------------------------------------------------*/

LevelMesh2d*
GascoigneMeshConstructor::LevelUpdate2d(GascoigneMesh* GM,
                                        const IndexSet& newquads,
                                        const IndexSet& oldquads) const
{
  LevelMesh2d* CM = new LevelMesh2d(HM);
  CM->BasicInit(newquads, oldquads);
  Construct2d(GM, CM);
  return CM;
}

/*---------------------------------------------------*/

LevelMesh3d*
GascoigneMeshConstructor::LevelUpdate3d(GascoigneMesh* GM,
                                        const IndexSet& newquads,
                                        const IndexSet& oldquads) const

{
  LevelMesh3d* CM = new LevelMesh3d(HM);
  CM->BasicInit(newquads, oldquads);
  Construct3d(GM, CM);
  return CM;
}

/*---------------------------------------------------*/

void
GascoigneMeshConstructor::Construct2d(GascoigneMesh* NNM,
                                      const LevelMesh2d* LM) const
{
  GascoigneMesh2d* NM = dynamic_cast<GascoigneMesh2d*>(NNM);

  assert(NM);

  IndexVector& nc = NM->GetCellVector();
  vector<Vertex2d>& nx = NM->GetVertexVector();
  IndexVector& mat = NM->GetMaterialVector();
  IndexVector& matpatch = NM->GetMaterialPatchVector();

  // zellen
  nc.reservesize(4 * LM->ncells());
  mat.reservesize(LM->ncells());
  for (IndexType i = 0; i < LM->ncells(); i++) {
    for (IndexType ii = 0; ii < 4; ii++)
      nc[4 * i + ii] = LM->vertex_of_cell(i, ii);
    mat[i] = LM->quad(i).material();
  }

  // Koordinaten

  nx.reserve(LM->nnodes());
  nx.resize(LM->nnodes());
  for (IndexType i = 0; i < LM->nnodes(); i++) {
    nx[i] = LM->vertex2d(i);
  }

  // PatchStructur

  PatchIndexHandler& PIH = NM->GetPatchIndexHandler();
  LM->ConstructIndOfPatch(PIH.GetIndex());
  PIH.GetDim() = 2;
  PIH.GetHasPatch() = 1;
  PatchToCell2d(PIH, LM);

  // Material Patch
  const nvector<IndexVector>& p2c = PIH.GetAllPatch2Cell();
  assert(p2c.size() == NM->npatches());
  matpatch.resize(p2c.size());
  for (IndexType ip = 0; ip < NM->npatches(); ++ip) {
    assert(p2c[ip].size() > 0);
    assert(p2c[ip][0] < mat.size());
    matpatch[ip] = mat[p2c[ip][0]];
  }
  // BoundaryIndices

  LM->InitBoundaryHandler(NNM->GetBoundaryIndexHandler(), PIH);

  // Hanging nodes

  HangingIndexHandler& H = NNM->GetHangingIndexHandler();
  LM->ConstructHangingStructureQuadratic(*(H.GetStructure()));
}

/*-----------------------------------------*/

void
GascoigneMeshConstructor::Construct3d(GascoigneMesh* NNM,
                                      const LevelMesh3d* LM) const
{
  GascoigneMesh3d* NM = dynamic_cast<GascoigneMesh3d*>(NNM);

  assert(NM);

  IndexVector& nc = NM->GetCellVector();
  vector<Vertex3d>& nx = NM->GetVertexVector();
  IndexVector& mat = NM->GetMaterialVector();
  IndexVector& matpatch = NM->GetMaterialPatchVector();

  IndexVector& mat_Vanka = NM->GetMaterialVankaVector();
  IndexVector& matpatch_Vanka = NM->GetMaterialVankaPatchVector();

  std::vector<std::array<Vertex3d, 3>>& basis_Vanka = NM->GetbasisVankaVector();
  std::vector<std::array<Vertex3d, 3>>& basis_Vanka_patch =
    NM->GetbasisVankaPatchVector();
  // zellen

  nc.reservesize(8 * LM->ncells());
  mat.reservesize(LM->ncells());
  mat_Vanka.reservesize(LM->ncells());
  basis_Vanka.resize(LM->ncells());

  for (IndexType i = 0; i < LM->ncells(); i++) {
    for (IndexType ii = 0; ii < 8; ii++)
      nc[8 * i + ii] = LM->vertex_of_cell(i, ii);
    mat[i] = LM->hex(i).material();
    mat_Vanka[i] = LM->hex(i).material_Vanka();
    basis_Vanka[i] = LM->hex(i).basis_Vanka();
  }

  // Koordinaten

  nx.reserve(LM->nnodes());
  nx.resize(LM->nnodes());
  for (IndexType i = 0; i < LM->nnodes(); i++) {
    nx[i] = LM->vertex3d(i);
  }

  // PatchStructur

  PatchIndexHandler& PIH = NM->GetPatchIndexHandler();

  LM->ConstructIndOfPatch(PIH.GetIndex());

  PIH.GetDim() = 3;
  PIH.GetHasPatch() = 1;

  PatchToCell3d(PIH, LM);

  // Material Patch
  const nvector<IndexVector>& p2c = PIH.GetAllPatch2Cell();
  assert(p2c.size() == NM->npatches());
  matpatch.resize(p2c.size());
  matpatch_Vanka.resize(p2c.size());
  basis_Vanka_patch.resize(p2c.size());

  for (IndexType ip = 0; ip < NM->npatches(); ++ip) {
    assert(p2c[ip].size() > 0);
    assert(p2c[ip][0] < mat.size());
    matpatch[ip] = mat[p2c[ip][0]];
    matpatch_Vanka[ip] = mat_Vanka[p2c[ip][0]];
    basis_Vanka_patch[ip] = basis_Vanka[p2c[ip][0]];
  }
  // BoundaryIndices

  LM->InitBoundaryHandler(NNM->GetBoundaryIndexHandler(), PIH);

  // Hanging nodes

  HangingIndexHandler& H = NNM->GetHangingIndexHandler();
  LM->ConstructHangingStructureQuadratic(*(H.GetStructure()),
                                         *(H.GetStructureFace()));
}

/*-----------------------------------------*/

void
GascoigneMeshConstructor::PatchToCell2d(PatchIndexHandler& PIH,
                                        const LevelMesh2d* LM) const
{
  IndexVector ci;
  LM->ConstructCellIndOfPatch(ci);

  IndexType np = ci.size();

  nvector<IndexVector>& patch2cell = PIH.GetAllPatch2Cell();
  patch2cell.resize(np);

  const HierarchicalMesh2d* HM2d = dynamic_cast<const HierarchicalMesh2d*>(HM);
  assert(HM2d);

  for (IndexType i = 0; i < np; ++i) {
    IndexType j = ci[i];
    const Quad& q = HM2d->quad(j);
    patch2cell[i].resize(4);
    for (IndexType p = 0; p < 4; ++p) {
      patch2cell[i][p] = LM->Quadg2l(q.child(p));
      assert(patch2cell[i][p] >= 0);
    }
  }
}

/*-----------------------------------------*/

void
GascoigneMeshConstructor::PatchToCell3d(PatchIndexHandler& PIH,
                                        const LevelMesh3d* LM) const
{
  IndexVector ci;
  LM->ConstructCellIndOfPatch(ci);

  IndexType np = ci.size();

  nvector<IndexVector>& patch2cell = PIH.GetAllPatch2Cell();
  patch2cell.resize(np);

  const HierarchicalMesh3d* HM3d = dynamic_cast<const HierarchicalMesh3d*>(HM);
  assert(HM3d);

  for (IndexType i = 0; i < np; ++i) {
    IndexType j = ci[i];
    const Hex& q = HM3d->hex(j);
    patch2cell[i].resize(8);
    for (IndexType p = 0; p < 8; ++p) {
      patch2cell[i][p] = LM->Hexg2l(q.child(p));
      assert(patch2cell[i][p] >= 0);
    }
  }
}
} // namespace Gascoigne

/*-----------------------------------------*/
