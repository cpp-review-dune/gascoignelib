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

#include "baseq12d.h"
#include "finiteelement.h"
#include "finiteelement.xx"
#include "transformation2d.h"

/*-----------------------------------------------------*/

namespace Gascoigne {
typedef Transformation2d<BaseQ12d> TQ1_2D;
template class FiniteElement<2, 1, Transformation2d<BaseQ12d>, BaseQ12d>;
} // namespace Gascoigne

/*-----------------------------------------------------*/
