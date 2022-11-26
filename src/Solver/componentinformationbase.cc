/**
 *
 * Copyright (C) 2005, 2011 by the Gascoigne 3D authors
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

#include "componentinformationbase.h"

#include "../Common/compose_name.h"
#include "../Interface/domainrighthandside.h"
#include "../Interface/problemdescriptorinterface.h"

namespace Gascoigne {
std::string
ComponentInformationBase::GetName() const
{
  return "ComponentInformationBase";
}

const IndexType
ComponentInformationBase::GetNScalars() const
{
  ProblemDescriptorInterface* ppdi = GetProblemDescriptorInterface();
  assert(ppdi);

  return ppdi->GetNcomp();
}
void
ComponentInformationBase::GetScalarName(IndexType i, std::string& s_name) const
{
  s_name = "u";
  compose_name_without_dot(s_name, i);
}
const IndexType
ComponentInformationBase::GetNVectors() const
{
  int ncomps = GetNcomp();
  if (ncomps <= 2)
    return 0;
  return 1;
}
void
ComponentInformationBase::GetVectorName(IndexType i, std::string& s_name) const
{
  s_name = "v";
}
void
ComponentInformationBase::GetVectorIndices(
  IndexType i,
  std::array<int, 3>& fa_vectorindices) const
{
  if (GetDimension() == 2) {
    fa_vectorindices[0] = 1;
    fa_vectorindices[1] = 2;
    fa_vectorindices[2] = -1;
  } else if (GetDimension() == 3) {
    fa_vectorindices[0] = 1;
    fa_vectorindices[1] = 2;
    fa_vectorindices[2] = 3;
  } else {
    std::cerr << __FILE__ << " :bad dimension=" << GetDimension() << "."
              << std::endl;
    abort();
  }
}

} // namespace Gascoigne
