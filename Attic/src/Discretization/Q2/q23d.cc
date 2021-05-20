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

#include "q23d.h"
#include "baseq23d.h"
#include "finiteelement.h"
#include "galerkinintegratorq2.h"
#include "gascoignemeshtransfer.h"
#include "hnstructureq23d.h"
#include "mginterpolatormatrix.h"
#include "mginterpolatornested.h"
#include "sparsestructure.h"
#include "transformation3d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne {
Q23d::Q23d()
  : Q2()
{
  HN = new HNStructureQ23d;
}

/* ----------------------------------------- */

Q23d::~Q23d()
{
  if (HN)
    delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

/* ----------------------------------------- */

void
Q23d::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
{
  const IntVector& vo2n = *GetMesh()->Vertexo2n();
  nvector<bool> habschon(GetMesh()->nnodes(), 0);

  assert(vo2n.size() == uold.n());

  for (int i = 0; i < vo2n.size(); i++) {
    int in = vo2n[i];

    if (in >= 0) {
      u.equ_node(in, 1., i, uold);
      habschon[in] = 1;
    }
  }
  nvector<std::array<int, 3>> nodes(12);
  nodes[0][0] = 1;
  nodes[0][1] = 0;
  nodes[0][2] = 2;
  nodes[1][0] = 3;
  nodes[1][1] = 0;
  nodes[1][2] = 6;
  nodes[2][0] = 5;
  nodes[2][1] = 2;
  nodes[2][2] = 8;
  nodes[3][0] = 7;
  nodes[3][1] = 6;
  nodes[3][2] = 8;
  nodes[4][0] = 1 + 18;
  nodes[4][1] = 0 + 18;
  nodes[4][2] = 2 + 18;
  nodes[5][0] = 3 + 18;
  nodes[5][1] = 0 + 18;
  nodes[5][2] = 6 + 18;
  nodes[6][0] = 5 + 18;
  nodes[6][1] = 2 + 18;
  nodes[6][2] = 8 + 18;
  nodes[7][0] = 7 + 18;
  nodes[7][1] = 6 + 18;
  nodes[7][2] = 8 + 18;
  nodes[8][0] = 9;
  nodes[8][1] = 0;
  nodes[8][2] = 18;
  nodes[9][0] = 11;
  nodes[9][1] = 2;
  nodes[9][2] = 20;
  nodes[10][0] = 15;
  nodes[10][1] = 6;
  nodes[10][2] = 24;
  nodes[11][0] = 17;
  nodes[11][1] = 8;
  nodes[11][2] = 26;

  nvector<std::array<int, 5>> w(6);
  w[0][0] = 4;
  w[0][1] = 0;
  w[0][2] = 2;
  w[0][3] = 6;
  w[0][4] = 8;
  w[1][0] = 12;
  w[1][1] = 0;
  w[1][2] = 18;
  w[1][3] = 6;
  w[1][4] = 24;
  w[2][0] = 14;
  w[2][1] = 2;
  w[2][2] = 8;
  w[2][3] = 20;
  w[2][4] = 26;
  w[3][0] = 16;
  w[3][1] = 6;
  w[3][2] = 8;
  w[3][3] = 24;
  w[3][4] = 26;
  w[4][0] = 10;
  w[4][1] = 0;
  w[4][2] = 2;
  w[4][3] = 18;
  w[4][4] = 20;
  w[5][0] = 22;
  w[5][1] = 18;
  w[5][2] = 20;
  w[5][3] = 24;
  w[5][4] = 26;

  for (int iq = 0; iq < GetMesh()->npatches(); ++iq) {
    IntVector vi = *GetMesh()->IndicesOfPatch(iq);

    for (int j = 0; j < nodes.size(); j++) {
      int v = vi[nodes[j][0]];
      int v1 = vi[nodes[j][1]];
      int v2 = vi[nodes[j][2]];
      assert(habschon[v1]);
      assert(habschon[v2]);
      if (habschon[v] == 0) {
        u.equ_node(v, 0.5, v1, uold);
        u.add_node(v, 0.5, v2, uold);
        habschon[v] = 1;
      }
    }
    for (int j = 0; j < w.size(); j++) {
      int v = vi[w[j][0]];
      int v1 = vi[w[j][1]];
      int v2 = vi[w[j][2]];
      int v3 = vi[w[j][3]];
      int v4 = vi[w[j][4]];
      assert(habschon[v1]);
      assert(habschon[v2]);
      assert(habschon[v3]);
      assert(habschon[v4]);
      if (habschon[v] == 0) {
        u.equ_node(v, 0.25, v1, uold);
        u.add_node(v, 0.25, v2, uold);
        u.add_node(v, 0.25, v3, uold);
        u.add_node(v, 0.25, v4, uold);
        habschon[v] = 1;
      }
    }
    int v = vi[13];
    if (habschon[v] == 0) {
      u.equ_node(v, 0.125, vi[0], uold);
      u.add_node(v, 0.125, vi[2], uold);
      u.add_node(v, 0.125, vi[6], uold);
      u.add_node(v, 0.125, vi[8], uold);
      u.add_node(v, 0.125, vi[18], uold);
      u.add_node(v, 0.125, vi[20], uold);
      u.add_node(v, 0.125, vi[24], uold);
      u.add_node(v, 0.125, vi[26], uold);
      habschon[v] = 1;
    }
  }
}

int
Q23d::GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const
{
  int iq;

  for (iq = 0; iq < GetMesh()->npatches(); ++iq) {
    bool found = true;
    const IntVector& IOP = GetMesh()->CoarseIndices(iq);

    for (int d = 0; d < 3; ++d) {
      double min = GetMesh()->vertex3d(IOP[0])[d];
      double max = min;
      for (int j = 1; j < 8; ++j) {
        double x = GetMesh()->vertex3d(IOP[j])[d];

        min = std::min(min, x);
        max = std::max(max, x);
      }
      if ((p0[d] < min) || (p0[d] > max)) {
        found = false;
        break;
      }
    }

    if (!found) {
      continue;
    }

    VertexTransformation(p0, p, iq);

    for (int d = 0; d < 3; ++d) {
      if ((p[d] < 0. - 1.e-12) || (p[d] > 1. + 1.e-12)) {
        found = false;
      }
    }

    if (found) {
      break;
    }
  }

  if (iq < GetMesh()->npatches()) {
    return iq;
  } else {
    return -1;
  }
}

/* ----------------------------------------- */

void
Q23d::VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T, iq);

  Transformation3d<BaseQ23d> Tr;
  Tr.init(T);

  Vertex3d res;

  p = 0.5;

  for (int niter = 1;; niter++) {
    Tr.point(p);

    res = p0;
    res.add(-1, Tr.x());

    if (res.norm() < 1.e-13) {
      break;
    }
    assert(niter < 10);

    Tr.DTI().mult_ad(p, res);
  }
}

/* ----------------------------------------- */

void
Q23d::BasicInit(const ParamFile* paramfile)
{
  if (!PatchDiscretization::GetIntegrator())
    PatchDiscretization::GetIntegratorPointer() = new GalerkinIntegratorQ2<3>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  if (!PatchDiscretization::GetFem()) {
    typedef Transformation3d<BaseQ23d> TransQ2;
    typedef FiniteElement<3, 2, TransQ2, BaseQ23d> FiniteElement;

    PatchDiscretization::GetFemPointer() = new FiniteElement;
  }
  assert(GetFem());

  PatchDiscretization::BasicInit(paramfile);
}

/* ----------------------------------------- */

void
Q23d::ConstructInterpolator(MgInterpolatorInterface* I,
                            const MeshTransferInterface* MT)
{
  MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);

  assert(IP);
  IP->BasicInit(MT);
  return;
}
} // namespace Gascoigne
