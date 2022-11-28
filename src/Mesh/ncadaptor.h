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

#ifndef __numerusclaususadaptor
#define __numerusclaususadaptor

#include <string>

#include "../Common/paramfile.h"
#include "../Interface/gascoigne.h"

/*-----------------------------------------------------------*/

namespace Gascoigne {
class NCAdaptor
{
protected:
  IndexType _n;
  double _p;
  double etasum;
  const DoubleVector& eta;

public:
  NCAdaptor(const ParamFile& paramfile, const DoubleVector& eta);
  void refine(IntVector& ref, IntVector& coarse) const;
};
} // namespace Gascoigne

#endif
