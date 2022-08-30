/**
 *
 * Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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

#include "gascoignemeshtransferconstructor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
GascoigneMeshTransferConstructor2d::GascoigneMeshTransferConstructor2d(
  const HierarchicalMesh2d* HM,
  GascoigneMeshTransfer* GMT,
  const LevelMesh2d* LMfine,
  const LevelMesh2d* LMcoarse)
{
  IntVector& c2f = GMT->GetC2f();
  map<IndexType, std::array<IndexType, 2>>& zweier = GMT->GetZweier();
  map<IndexType, std::array<IndexType, 4>>& vierer = GMT->GetVierer();
  map<IndexType, IndexType>& CellEiner = GMT->GetCellEiner();
  map<IndexType, std::array<IndexType, 4>>& CellVierer = GMT->GetCellVierer();

  const QuadLawAndOrder& LaO = HM->QuadLawOrder();

  CellEiner.clear();
  CellVierer.clear();
  // Zellen
  for (IndexType i = 0; i < LMcoarse->ncells(); i++) {
    IndexType igq = LMcoarse->Quadl2g(i);
    if (LMfine->Quadg2l().find(igq) != LMfine->Quadg2l().end()) {
      CellEiner[i] = LMfine->Quadg2l(igq);
    } else {
      // kinder suchen
      std::array<IndexType, 4> n4;
      for (IndexType ii = 0; ii < 4; ii++) {
        IndexType ic = HM->quad(igq).child(ii);
        n4[ii] = LMfine->Quadg2l(ic);
      }
      CellVierer[i] = n4;
    }
  }

  // 1er
  c2f.reservesize(LMcoarse->nnodes());
  for (IndexType iL = 0; iL < LMcoarse->nnodes(); iL++) {
    IndexType ig = LMcoarse->Vertexl2g(iL);
    IndexType il = LMfine->Vertexg2l(ig);

    c2f[iL] = il;
  }

  // 2er 4er 8er (!)
  for (IndexType i = 0; i < LMcoarse->ncells(); i++) {
    IndexType igq = LMcoarse->Quadl2g(i);
    if (LMfine->Quadg2l().find(igq) != LMfine->Quadg2l().end())
      continue;
    const Quad& q = HM->quad(igq);

    // verfeinertes quad --> vierer

    std::array<IndexType, 4> n4;
    IndexType igm = LaO.middle_vertex(q);
    IndexType ilm = LMfine->Vertexg2l(igm);
    for (IndexType ii = 0; ii < 4; ii++) {
      IndexType ig = q.vertex(ii);
      IndexType iL = LMcoarse->Vertexg2l(ig);
      n4[ii] = iL;
    }
    vierer[ilm] = n4;

    // edges
    for (IndexType ie = 0; ie < 4; ie++) {
      IndexType ige = LaO.edge_vertex(q, ie);
      IndexType ile = LMfine->Vertexg2l(ige);
      if (LMcoarse->Vertexg2lCheck(ige) != -2)
        continue;

      std::array<IndexType, 2> f;
      LaO.globalvertices_of_edge(q, f, ie);

      std::array<IndexType, 2> n2;
      for (IndexType ii = 0; ii < 2; ii++) {
        IndexType iL = LMcoarse->Vertexg2l(f[ii]);
        n2[ii] = iL;
      }
      zweier[ile] = n2;
    }
  }
  //   ofstream file("MGI",ios::app);
  //   cout << "Matrix:\n" << M << endl;
}

/*-----------------------------------------*/

GascoigneMeshTransferConstructor3d::GascoigneMeshTransferConstructor3d(
  const HierarchicalMesh3d* HM,
  GascoigneMeshTransfer* GMT,
  const LevelMesh3d* LMfine,
  const LevelMesh3d* LMcoarse)
{
  //   cerr << "GascoigneMeshTransferConstructor::Construct3d()\n";
  //   cerr << "noch keine konstanten!\n";
  //   abort();

  IntVector& c2f = GMT->GetC2f();
  map<IndexType, std::array<IndexType, 2>>& zweier = GMT->GetZweier();
  map<IndexType, std::array<IndexType, 4>>& vierer = GMT->GetVierer();
  map<IndexType, std::array<IndexType, 8>>& achter = GMT->GetAchter();

  const HexLawAndOrder& LaO = HM->HexLawOrder();

  // 1er
  c2f.reservesize(LMcoarse->nnodes());
  for (IndexType iL = 0; iL < LMcoarse->nnodes(); iL++) {
    IndexType ig = LMcoarse->Vertexl2g(iL);
    IndexType il = LMfine->Vertexg2l(ig);
    c2f[iL] = il;
  }
  // 2er 4er 8er (!)
  for (IndexType i = 0; i < LMcoarse->ncells(); i++) {
    IndexType igq = LMcoarse->Hexl2g(i);
    if (LMfine->Hexg2l().find(igq) != LMfine->Hexg2l().end())
      continue;
    const Hex& q = HM->hex(igq);

    // verfeinertes hex

    std::array<IndexType, 8> n8;
    IndexType igm = LaO.middle_vertex(q);
    IndexType ilm = LMfine->Vertexg2l(igm);
    for (IndexType ii = 0; ii < 8; ii++) {
      IndexType ig = q.vertex(ii);
      // 8er
      n8[ii] = LMcoarse->Vertexg2l(ig);
    }
    achter[ilm] = n8;
    // faces
    for (IndexType ie = 0; ie < 6; ie++) {
      IndexType ige = LaO.face_vertex(q, ie);
      IndexType ile = LMfine->Vertexg2l(ige);
      if (LMcoarse->Vertexg2lCheck(ige) != -2)
        continue;
      std::array<IndexType, 4> f, n4;
      LaO.globalvertices_of_face(q, f, ie);
      for (IndexType ii = 0; ii < 4; ii++) {
        // 4er
        n4[ii] = LMcoarse->Vertexg2l(f[ii]);
      }
      vierer[ile] = n4;
    }
    // edges
    for (IndexType ie = 0; ie < 12; ie++) {
      IndexType ige = LaO.edge_vertex(q, ie);
      IndexType ile = LMfine->Vertexg2l(ige);
      if (LMcoarse->Vertexg2lCheck(ige) != -2)
        continue;
      std::array<IndexType, 2> f, n2;
      LaO.globalvertices_of_edge(q, f, ie);

      n2[0] = LMcoarse->Vertexg2l(f[0]);
      n2[1] = LMcoarse->Vertexg2l(f[1]);

      zweier[ile] = n2;
    }
  }
}
} // namespace Gascoigne
