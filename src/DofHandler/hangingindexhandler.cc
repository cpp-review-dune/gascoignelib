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

#include "hangingindexhandler.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
HangingIndexHandler::HangingIndexHandler()
{
  hnq2.clear();
  hnq2face.clear();
}

/*-----------------------------------------*/

void
HangingIndexHandler::CopyLevel2Nibble(const HangingIndexHandler& Lhih,
                                      IndexVector& Vg2l)
{
  hnq2.clear();
  hnq2face.clear();

  map<IndexType, std::array<IndexType, 3>>::const_iterator it3 =
    Lhih.GetStructure()->begin();
  map<IndexType, std::array<IndexType, 9>>::const_iterator it9 =
    Lhih.GetStructureFace()->begin();
  map<IndexType, std::array<IndexType, 3>>::const_iterator end3 =
    Lhih.GetStructure()->end();
  map<IndexType, std::array<IndexType, 9>>::const_iterator end9 =
    Lhih.GetStructureFace()->end();

  for (; it3 != end3; ++it3) {
    IndexType gf = it3->first;
    IndexType lf = Vg2l[gf];
    if (lf < 0)
      continue;
    std::array<IndexType, 3> tmp;
    IndexType gut = 3;
    for (IndexType i = 0; i < 3; ++i) {
      tmp[i] = Vg2l[it3->second[i]];
      if (tmp[i] < 0)
        --gut;
      if (i < 2)
        assert(tmp[i] >= 0);
    }
    assert(gut == 3);
    hnq2[lf] = tmp;
  }
  for (; it9 != end9; ++it9) {
    IndexType gf = it9->first;
    IndexType lf = Vg2l[gf];
    if (lf < 0)
      continue;
    std::array<IndexType, 9> tmp;
    IndexType gut = 9;
    for (IndexType i = 0; i < 9; ++i) {
      tmp[i] = Vg2l[it9->second[i]];
      if (tmp[i] < 0)
        --gut;
    }
    assert(gut == 9);
    hnq2face[lf] = tmp;
  }
}
} // namespace Gascoigne
