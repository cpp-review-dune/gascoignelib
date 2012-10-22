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


#include  "nodematrix.h"

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/*-------------------------------------------------------*/

namespace Gascoigne
{
template class NodeMatrix<1 ,float>;
template class NodeMatrix<2 ,float>;
template class NodeMatrix<3 ,float>;
template class NodeMatrix<4 ,float>;
template class NodeMatrix<5 ,float>;
template class NodeMatrix<6 ,float>;
template class NodeMatrix<7 ,float>;
template class NodeMatrix<8 ,float>;
template class NodeMatrix<9 ,float>;
template class NodeMatrix<16,float>;
template class NodeMatrix<20,float>;
template class NodeMatrix<25,float>;
}
