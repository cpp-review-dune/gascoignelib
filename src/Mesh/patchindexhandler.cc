/**
 *
 * Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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

#include "patchindexhandler.h"

/*-----------------------------------------*/

namespace Gascoigne {

void
PatchIndexHandler::ConstructQ4PatchByQ2Patch(
  const nvector<IntVector>& coarseq2patches)
{
  assert(haspatch);
  if (dim == 2) {
    assert(coarseq2patches.size() * 4 == indexofpatch.size());

    hasq4patch = true;

    // construct temporary map on fine mesh: from (n1,n2) -> q2patch
    // (n1,n2) are the two opposite nodes of the q2 patch.
    std::map<std::array<int, 2>, int> fine_diag2q2patch;
    for (int q = 0; q < indexofpatch.size(); ++q) {
      const IntVector& pi = indexofpatch[q];
      assert(pi.size() == 9); // only 2d
      fine_diag2q2patch[{ pi[0], pi[8] }] = q;
    }

    indexofq4patch.clear();
    indexofq4patch.resize(coarseq2patches.size(), IntVector(25, 0.));

    // run over all coarse q2 patches and enter indices and cells

    for (int q = 0; q < coarseq2patches.size(); ++q) {
      const IntVector& cpi = coarseq2patches[q];

      assert(cpi.size() == 9);

      for (int iy = 0; iy < 2; ++iy) {
        for (int ix = 0; ix < 2; ++ix) {
          int n0 = cpi[3 * iy + ix]; // lower/left node of patch
          int n1 = cpi[3 * iy + ix + 4];

          assert(fine_diag2q2patch.find({ n0, n1 }) != fine_diag2q2patch.end());
          const IntVector& fpi = indexofpatch[fine_diag2q2patch[{ n0, n1 }]];

          int F0 = 10 * iy + 2 * ix; // first index in fine patch
          for (int cy = 0; cy < 3; ++cy)
            for (int cx = 0; cx < 3; ++cx)
              indexofq4patch[q][F0 + 5 * cy + cx] = fpi[3 * cy + cx];
        }
      }
    }
  } else if (dim == 3) {
    assert(coarseq2patches.size() * 8 == indexofpatch.size());

    hasq4patch = true;

    // construct temporary map on fine mesh: from (n1,n2) -> q2patch
    // (n1,n2) are the two opposite nodes of the q2 patch.
    std::map<std::array<int, 2>, int> fine_diag2q2patch;
    for (int q = 0; q < indexofpatch.size(); ++q) {
      const IntVector& pi = indexofpatch[q];
      assert(pi.size() == 27); // only 2d
      fine_diag2q2patch[{ pi[0], pi[26] }] = q;
    }

    indexofq4patch.clear();
    indexofq4patch.resize(coarseq2patches.size(), IntVector(125, 0.));

    // run over all coarse q2 patches and enter indices and cells

    for (int q = 0; q < coarseq2patches.size(); ++q) {
      const IntVector& cpi = coarseq2patches[q];

      assert(cpi.size() == 27);
      for (int iz = 0; iz < 2; ++iz) {
        for (int iy = 0; iy < 2; ++iy) {
          for (int ix = 0; ix < 2; ++ix) {
            int n0 = cpi[9 * iz + 3 * iy + ix]; // lower/left node of patch
            int n1 = cpi[9 * iz + 3 * iy + ix + 13];
            assert(fine_diag2q2patch.find({ n0, n1 }) !=
                   fine_diag2q2patch.end());
            const IntVector& fpi = indexofpatch[fine_diag2q2patch[{ n0, n1 }]];

            int F0 = 50 * iz + 10 * iy + 2 * ix; // first index in fine patch
            for (int cz = 0; cz < 3; ++cz)
              for (int cy = 0; cy < 3; ++cy)
                for (int cx = 0; cx < 3; ++cx)
                  indexofq4patch[q][F0 + 25 * cz + 5 * cy + cx] =
                    fpi[9 * cz + 3 * cy + cx];
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////// TEMP 06. October 2021 can
/// all be removed

IntVector
PatchIndexHandler::Q2IndicesOfQ4Patch(int i) const
{
  assert(hasq4patch);

  IntVector indices;
  if (dim == 2) {
    indices.resize(9);
    indices[0] = indexofq4patch[i][0];
    indices[1] = indexofq4patch[i][2];
    indices[2] = indexofq4patch[i][4];
    indices[3] = indexofq4patch[i][10];
    indices[4] = indexofq4patch[i][12];
    indices[5] = indexofq4patch[i][14];
    indices[6] = indexofq4patch[i][20];
    indices[7] = indexofq4patch[i][22];
    indices[8] = indexofq4patch[i][24];
  } else {
    indices.resize(0);
    for (int z = 0; z < 3; z++) {
      indices.push_back(indexofq4patch[i][50 * z + 0]);
      indices.push_back(indexofq4patch[i][50 * z + 2]);
      indices.push_back(indexofq4patch[i][50 * z + 4]);
      indices.push_back(indexofq4patch[i][50 * z + 10]);
      indices.push_back(indexofq4patch[i][50 * z + 12]);
      indices.push_back(indexofq4patch[i][50 * z + 14]);
      indices.push_back(indexofq4patch[i][50 * z + 20]);
      indices.push_back(indexofq4patch[i][50 * z + 22]);
      indices.push_back(indexofq4patch[i][50 * z + 24]);
    }
  }
  return indices;
}

/*-----------------------------------------*/

IntVector
PatchIndexHandler::CoarseIndices(int iq) const
{
  IntVector indices;

  assert(iq >= 0);
  assert(iq < indexofpatch.size());

  if (dim == 2) {
    indices.resize(4);

    indices[0] = indexofpatch[iq][0];
    indices[1] = indexofpatch[iq][2];
    indices[2] = indexofpatch[iq][6];
    indices[3] = indexofpatch[iq][8];
  } else {
    indices.resize(8);

    indices[0] = indexofpatch[iq][0];
    indices[1] = indexofpatch[iq][2];
    indices[2] = indexofpatch[iq][6];
    indices[3] = indexofpatch[iq][8];
    indices[4] = indexofpatch[iq][18];
    indices[5] = indexofpatch[iq][20];
    indices[6] = indexofpatch[iq][24];
    indices[7] = indexofpatch[iq][26];
  }
  return indices;
}

/*-----------------------------------------*/

IntVector
PatchIndexHandler::CoarseIndicesQ4(int iq) const
{
  assert(hasq4patch);

  IntVector indices;

  assert(iq >= 0);
  assert(iq < indexofq4patch.size());

  if (dim == 2) {
    indices.resize(4);

    indices[0] = indexofq4patch[iq][0];
    indices[1] = indexofq4patch[iq][4];
    indices[2] = indexofq4patch[iq][20];
    indices[3] = indexofq4patch[iq][24];
  } else {
    indices.resize(8);

    indices[0] = indexofq4patch[iq][0];
    indices[1] = indexofq4patch[iq][4];
    indices[2] = indexofq4patch[iq][20];
    indices[3] = indexofq4patch[iq][24];
    indices[4] = indexofq4patch[iq][100];
    indices[5] = indexofq4patch[iq][104];
    indices[6] = indexofq4patch[iq][120];
    indices[7] = indexofq4patch[iq][124];
  }
  return indices;
}
} // namespace Gascoigne
