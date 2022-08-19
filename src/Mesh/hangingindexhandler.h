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

#ifndef __HangingIndexHandler_h
#define __HangingIndexHandler_h

#include "gascoigne.h"
#include <array>
#include <map>

/*-----------------------------------------*/

namespace Gascoigne {
class HangingIndexHandler
{
protected:
  // typedef std::array<IndexType, 2> IntVector2;
  // typedef std::array<IndexType, 4> IntVector4;

  std::map<IndexType, std::array<IndexType, 3>> hnq2;

  // // Average
  // U[it->first] = 0.5 * (U[it->second[0]] + U[it->second[1]]);
  // // Zero
  // U[it->first] = 0.
  // // Distribute
  // U[it->second[0] ] += 0.5 * U[it->first]
  // U[it->second[1] ] += 0.5 * U[it->first]

  // -S0-*--*
  //  *  |  |
  //  fa-*--*
  //  *  |  |
  // -S1-*- *
  //  *  |  |
  //  fb-*--*
  //  *  |  |
  // -S2 *- *

  //  fa -> (S0,S1,S2)
  //  fb -> (S2,S1,S0)

  // fuer 3D
  //
  std::map<IndexType, std::array<IndexType, 9>> hnq2face;
  // // Average
  // U[it->first] = 0.25 * (U[it->second[0]] + U[it->second[1]] +
  // U[it->second[3]] + U[it->second[4]]);
  // // Zero
  // U[it->first] = 0.
  // // Distribute
  // U[it->second[0] ] += 0.25 * U[it->first]
  // U[it->second[1] ] += 0.25 * U[it->first]
  // U[it->second[3] ] += 0.25 * U[it->first]
  // U[it->second[4] ] += 0.25 * U[it->first]

  std::map<IndexType, std::array<IndexType, 6>> hnq4;
  std::map<IndexType, std::array<IndexType, 26>> hnq4face;

public:
  HangingIndexHandler();

  void Equal(const std::map<IndexType, std::array<IndexType, 3>>& h2,
             const std::map<IndexType, std::array<IndexType, 9>>& h2f)
  {
    hnq2 = h2;
    hnq2face = h2f;
  }

  void CopyLevel2Nibble(const HangingIndexHandler& Lhih, IntVector& Vg2l);

  // zugriff

  const std::map<IndexType, std::array<IndexType, 3>>* GetStructure() const
  {
    return &hnq2;
  }
  const std::map<IndexType, std::array<IndexType, 9>>* GetStructureFace() const
  {
    return &hnq2face;
  }
  std::map<IndexType, std::array<IndexType, 3>>* GetStructure()
  {
    return &hnq2;
  }
  std::map<IndexType, std::array<IndexType, 9>>* GetStructureFace()
  {
    return &hnq2face;
  }

  const std::map<IndexType, std::array<IndexType, 6>>* GetQ4Structure() const
  {
    return &hnq4;
  }
  const std::map<IndexType, std::array<IndexType, 26>>* GetQ4StructureFace()
    const
  {
    return &hnq4face;
  }
  std::map<IndexType, std::array<IndexType, 6>>* GetQ4Structure()
  {
    return &hnq4;
  }
  std::map<IndexType, std::array<IndexType, 26>>* GetQ4StructureFace()
  {
    return &hnq4face;
  }
};
} // namespace Gascoigne

#endif
