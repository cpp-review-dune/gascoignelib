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

#ifndef __GascoigneMeshTransfer_h
#define __GascoigneMeshTransfer_h

#include "gascoigne.h"
#include "meshtransferinterface.h"
#include <array>
#include <map>

/*-----------------------------------------*/

namespace Gascoigne {
class GascoigneMeshTransfer : public MeshTransferInterface
{
protected:
  std::map<IndexType, std::array<IndexType, 2>> zweier;
  std::map<IndexType, std::array<IndexType, 4>> vierer;
  std::map<IndexType, std::array<IndexType, 8>> achter;

  IntVector c2f;
  std::map<IndexType, IndexType> CellEiner;
  std::map<IndexType, std::array<IndexType, 4>> CellVierer;
  std::map<IndexType, std::array<IndexType, 8>> CellAchter;

public:
  const std::map<IndexType, std::array<IndexType, 2>>& GetZweier() const
  {
    return zweier;
  }
  const std::map<IndexType, std::array<IndexType, 4>>& GetVierer() const
  {
    return vierer;
  }
  const std::map<IndexType, std::array<IndexType, 8>>& GetAchter() const
  {
    return achter;
  }
  const IntVector& GetC2f() const { return c2f; }

  std::map<IndexType, std::array<IndexType, 2>>& GetZweier() { return zweier; }
  std::map<IndexType, std::array<IndexType, 4>>& GetVierer() { return vierer; }
  std::map<IndexType, std::array<IndexType, 8>>& GetAchter() { return achter; }
  IntVector& GetC2f() { return c2f; }

  const std::map<IndexType, IndexType>& GetCellEiner() const
  {
    return CellEiner;
  }
  const std::map<IndexType, std::array<IndexType, 4>>& GetCellVierer() const
  {
    return CellVierer;
  }
  const std::map<IndexType, std::array<IndexType, 8>>& GetCellAchter() const
  {
    return CellAchter;
  }

  std::map<IndexType, IndexType>& GetCellEiner() { return CellEiner; }
  std::map<IndexType, std::array<IndexType, 4>>& GetCellVierer()
  {
    return CellVierer;
  }
  std::map<IndexType, std::array<IndexType, 8>>& GetCellAchter()
  {
    return CellAchter;
  }

  GascoigneMeshTransfer();
};
} // namespace Gascoigne

#endif
