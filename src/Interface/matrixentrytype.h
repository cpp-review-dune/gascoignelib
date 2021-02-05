/*----------------------------   matrixentrytype.h ---------------------------*/
/*      $Id:$                 */
#ifndef __matrixentrytype_H
#define __matrixentrytype_H
/*----------------------------   matrixentrytype.h ---------------------------*/

/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#ifdef __MATRIX_DOUBLE_PRECISION__
#define MatrixEntryType double
#else
#ifdef __MATRIX_SINGLE_PRECISION__
#define MatrixEntryType float
#else
#include <stophere> // either define __MATRIX_DOUBLE_PRECISION__ or __MATRIX_SINGLE_PRECISION
#endif
#endif

/*----------------------------   matrixentrytype.h ---------------------------*/
/* end of #ifndef __matrixentrytype_H */
#endif
/*----------------------------   matrixentrytype.h ---------------------------*/
