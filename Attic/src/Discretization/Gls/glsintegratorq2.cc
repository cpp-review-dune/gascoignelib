/**
 *
 * Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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

#include "glsintegratorq2.h"

using namespace std;

namespace Gascoigne {
/* ----------------------------------------- */

template <int DIM> GlsIntegratorQ2<DIM>::GlsIntegratorQ2() {
  if (DIM == 2)
    GlsIntegrator<DIM>::IF = new QuadGauss9;
  else
    GlsIntegrator<DIM>::IF = new HexGauss27;
  assert(GlsIntegrator<DIM>::IF);
}

/*-----------------------------------------------------------*/

template class GlsIntegratorQ2<2>;
template class GlsIntegratorQ2<3>;
} // namespace Gascoigne
