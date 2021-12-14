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

#ifndef __GalerkinIntegratorQ2_h
#define __GalerkinIntegratorQ2_h

#include "galerkinintegrator.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinIntegratorQ2

////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinIntegratorQ2 : virtual public GalerkinIntegrator<DIM>
{
public:
  //
  ////  Con(De)structor
  //
  GalerkinIntegratorQ2();

  ~GalerkinIntegratorQ2();
  std::string GetName() const { return "GalerkinIntegratorQ2"; }

  void BasicInit();
};
} // namespace Gascoigne

#endif
