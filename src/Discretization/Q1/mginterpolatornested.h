/**
 *
 * Copyright (C) 2004, 2007 by the Gascoigne 3D authors
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

#ifndef __MgInterpolatorNested_h
#define __MgInterpolatorNested_h

#include "gascoigne.h"
#include "meshtransferinterface.h"
#include "mginterpolatorinterface.h"
#include <map>

/*-----------------------------------------*/

namespace Gascoigne {
class MgInterpolatorNested : public virtual MgInterpolatorInterface
{
private:
  std::map<int, std::array<int, 2>> zweier;
  std::map<int, std::array<int, 4>> vierer;
  std::map<int, std::array<int, 8>> achter;

  IntVector c2f;

public:
  std::map<int, std::array<int, 2>>& GetZweier() { return zweier; }
  std::map<int, std::array<int, 4>>& GetVierer() { return vierer; }
  std::map<int, std::array<int, 8>>& GetAchter() { return achter; }

  IntVector& GetC2F() { return c2f; }
  const IntVector& GetC2F() const { return c2f; }

  MgInterpolatorNested()
    : MgInterpolatorInterface()
  {}

  virtual void BasicInit(const MeshTransferInterface* MT);

  void restrict_zero(GlobalVector&, const GlobalVector&) const;
  void prolongate_add(GlobalVector&, const GlobalVector&) const;
  void SolutionTransfer(GlobalVector&, const GlobalVector&) const;
  void Pi(GlobalVector& u) const;
};
} // namespace Gascoigne

#endif
