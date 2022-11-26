/**
 *
 * Copyright (C) 2005 by the Gascoigne 3D authors
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

#ifndef __ComponentInformationBase_h
#define __ComponentInformationBase_h

#include "../Common/filescanner.h"
#include "../Common/gostream.h"
#include "../Common/stringutil.h"
#include "../Interface/componentinformation.h"

namespace Gascoigne {

/////////////////////////////////////////////
///
///@brief
///  ... comments ComponentInformationBase

///
///
/////////////////////////////////////////////

class ComponentInformationBase : public ComponentInformation
{
private:
protected:
public:
  ComponentInformationBase()
    : ComponentInformation()
  {}
  virtual ~ComponentInformationBase() {}

  virtual void BasicInit(const ParamFile* pf) {}

  virtual std::string GetName() const;

  virtual const IndexType GetNScalars() const;
  virtual void GetScalarName(IndexType i, std::string& s_name) const;
  virtual const IndexType GetNVectors() const;
  virtual void GetVectorName(IndexType i, std::string& s_name) const;
  virtual void GetVectorIndices(IndexType i,
                                std::array<int, 3>& fa_vectorindices) const;
};
} // namespace Gascoigne

#endif // __ComponentInformationBase_h
