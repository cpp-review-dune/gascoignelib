/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2010, 2018 by the Gascoigne 3D authors
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


#include "elementintegrator.xx"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{


  template class ElementIntegratorQ12d;
  template class ElementIntegratorQ22d;
  template class ElementIntegratorQ42d;
  
  template class ElementIntegratorQ13d;
  template class ElementIntegratorQ23d;

  ////////////////////////////////////////////////// required for LPS  
  template class ElementIntegrator<2, PatchFormula2d<4,QuadGauss4>, PatchFormula2d<9,QuadGauss9>,  PatchFormula1d<2,LineGauss2>, PatchFormula2d<4,QuadGauss4>>;
  template class ElementIntegrator<3, PatchFormula3d<8,HexGauss8>,  PatchFormula3d<27,HexGauss27>, PatchFormula2d<4,QuadGauss4>, PatchFormula3d<8,HexGauss8>>;

  
  /* ----------------------------------------- */

} // namespace Gascoigne
