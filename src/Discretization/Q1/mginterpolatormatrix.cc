/**
 *
 * Copyright (C) 2004 by the Gascoigne 3D authors
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

#include "mginterpolatormatrix.h"
#include "gascoignemeshtransfer.h"
#include "sparsestructure.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {

void
MgInterpolatorMatrix::BasicInit(const MeshTransferInterface* MT)
{
  assert(MT);
  const GascoigneMeshTransfer* GT =
    dynamic_cast<const GascoigneMeshTransfer*>(MT);
  assert(GT);

  const std::map<int, std::array<int, 2>>& I2 = GT->GetZweier();
  const std::map<int, std::array<int, 4>>& I4 = GT->GetVierer();
  const std::map<int, std::array<int, 8>>& I8 = GT->GetAchter();
  const IntVector& IC = GT->GetC2f();

  int ncoarse = IC.size();
  int nfine =
    -1; // number of fine nodes must be reconstructed from Mesh Transfer
  for (IndexType i = 0; i < IC.size(); ++i)
    nfine = std::max(IC[i], nfine);
  for (const auto i1 : I2)
    nfine = std::max(i1.first, nfine);
  for (const auto i1 : I4)
    nfine = std::max(i1.first, nfine);
  for (const auto i1 : I8)
    nfine = std::max(i1.first, nfine);
  nfine++;

  SparseStructure Sfine, Scoarse;
  Sfine.build_begin(nfine);
  Scoarse.build_begin(ncoarse);
  for (int i = 0; i < ncoarse; ++i) {
    Sfine.build_add(IC[i], i);
    Scoarse.build_add(i, IC[i]);
  }
  for (const auto i1 : I2)
    for (const auto i2 : i1.second) {
      Sfine.build_add(i1.first, i2);
      Scoarse.build_add(i2, i1.first);
    }

  for (const auto i1 : I4)
    for (const auto i2 : i1.second) {
      Sfine.build_add(i1.first, i2);
      Scoarse.build_add(i2, i1.first);
    }

  for (const auto i1 : I8)
    for (const auto i2 : i1.second) {
      Sfine.build_add(i1.first, i2);
      Scoarse.build_add(i2, i1.first);
    }

  Sfine.build_end();
  Scoarse.build_end();
  STfine.memory(&Sfine);
  STcoarse.memory(&Scoarse);

  valfine.clear();
  valcoarse.clear();
  valfine.resize(STfine.nentries(), 0.0);
  valcoarse.resize(STcoarse.nentries(), 0.0);

  for (IndexType row = 0; row < STfine.n(); ++row) {
    for (IndexType pos = STfine.start(row); pos < STfine.stop(row); ++pos) {
      if (I2.find(row) != I2.end())
        valfine[pos] = 0.5;
      else if (I4.find(row) != I4.end())
        valfine[pos] = 0.25;
      else if (I8.find(row) != I8.end())
        valfine[pos] = 0.125;
      else
        valfine[pos] = 1.0; // in I2C - List
    }
  }

  for (IndexType row = 0; row < STcoarse.n(); ++row) {
    for (IndexType pos = STcoarse.start(row); pos < STcoarse.stop(row); ++pos) {
      IndexType col = STcoarse.col(pos);
      if (I2.find(col) != I2.end())
        valcoarse[pos] = 0.5;
      else if (I4.find(col) != I4.end())
        valcoarse[pos] = 0.25;
      else if (I8.find(col) != I8.end())
        valcoarse[pos] = 0.125;
      else
        valcoarse[pos] = 1.0; // in I2C - List
    }
  }
}

void
MgInterpolatorMatrix::restrict_zero(GlobalVector& uL,
                                    const GlobalVector& ul) const
{
#pragma omp parallel for
  for (int i = 0; i < STcoarse.n(); i++) {
    uL.zero_node(i);
    for (int pos = STcoarse.start(i); pos < STcoarse.stop(i); pos++)
      uL.add_node(i, valcoarse[pos], STcoarse.col(pos), ul);
  }
}

/*-----------------------------------------*/

void
MgInterpolatorMatrix::prolongate_add(GlobalVector& ul,
                                     const GlobalVector& uL) const
{
#pragma omp parallel for
  for (int i = 0; i < STfine.n(); i++) {
    for (int pos = STfine.start(i); pos < STfine.stop(i); pos++) {
      ul.add_node(i, valfine[pos], STfine.col(pos), uL);
    }
  }
}

/*-----------------------------------------*/

void
MgInterpolatorMatrix::SolutionTransfer(GlobalVector& uL,
                                       const GlobalVector& ul) const
{
#pragma omp parallel for
  for (int i = 0; i < STcoarse.n(); i++) {
    uL.zero_node(i);
    for (int pos = STcoarse.start(i); pos < STcoarse.stop(i); pos++) {
      if (valcoarse[pos] == 1.)
        uL.add_node(i, 1.0, STcoarse.col(pos), ul);
    }
  }
}
} // namespace Gascoigne
