/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#include "finiteelement.h"
#include "../Q1/finiteelement.xx"
#include "baseq1.h"
#include "baseq2.h"
#include "baseq42d.h"
#include "baseq43d.h"
#include "baseq2patch.h"
#include "baseq23dpatch.h"
#include "transformation.h"

namespace Gascoigne
{

/**********************************************************/

  typedef Gascoigne::Transformation<2, BaseQ1<2> > TQ1_2D;
  typedef Gascoigne::Transformation<2, BaseQ2<2> > TQ2_2D;
  typedef Gascoigne::Transformation<2, BaseQ42d> TQ4_2D;

  template class FiniteElement<2,1,TQ1_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ2_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ4_2D,BaseQ42d>;
  template class FiniteElement<2,1,TQ1_2D,BaseQ22dPatch>;
  template class FiniteElement<2,1,TQ2_2D,BaseQ22dPatch>;

/**********************************************************/

  typedef Gascoigne::Transformation<3, BaseQ1<3> > TQ1_3D;
  typedef Gascoigne::Transformation<3, BaseQ2<3> > TQ2_3D;
  typedef Gascoigne::Transformation<3, BaseQ43d> TQ4_3D;

  template class FiniteElement<3,2,TQ1_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ2_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ4_3D,BaseQ43d>;
  template class FiniteElement<3,2,TQ1_3D,BaseQ23dPatch>;
  template class FiniteElement<3,2,TQ2_3D,BaseQ23dPatch>;

/**********************************************************/

}

