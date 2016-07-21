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


#include  "transformation.h"
#include  "baseq23dwithsecond.h"
#include  "../Q1/finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond3d.xx"

/*-----------------------------------------------------*/
namespace Gascoigne
{
  template class FiniteElement          <3, 2, Gascoigne::Transformation<3, BaseQ23dWithSecond>, BaseQ23dWithSecond>;
  template class FiniteElementWithSecond<3, 2, Gascoigne::Transformation<3, BaseQ23dWithSecond>, BaseQ23dWithSecond>;
}
/*-----------------------------------------------------*/
