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


#include  "finiteelement.h"
#include  "finiteelement.xx"
#include  "transformation2d.h"
#include  "baseq1.h"

/*-----------------------------------------------------*/

namespace Gascoigne
{
  typedef Transformation2d<BaseQ1<2> >  TQ1_2D;

  template class FiniteElement<2,1,TQ1_2D,BaseQ1<2> >;
}

/*-----------------------------------------------------*/
