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

#ifndef __Gascoigne_h
#define __Gascoigne_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Gascoigne

////
////
/////////////////////////////////////////////

#include <map>
#include <set>
#include <string>

#include "derivativevector.h"
#include "nmatrix.h"

namespace Gascoigne {

#ifdef __MATRIX_SINGLE_PRECISION__
typedef float MatrixEntryType;
#else
typedef double MatrixEntryType;
#endif

typedef size_t IndexType;                  // data type for global indices
typedef unsigned short int ShortIndexType;                  // data type for local indices
  

typedef nvector<double> GlobalParameterVector;
typedef nvector<double> LocalParameterVector;
typedef std::map<std::string, const GlobalParameterVector *>
    GlobalParameterData;
typedef std::map<std::string, LocalParameterVector> LocalParameterData;

typedef nvector<int> IntVector;
typedef nvector<IndexType> IndexVector;
typedef nvector<double> DoubleVector;
typedef nmatrix<double> DoubleMatrix;
typedef std::set<int> IntSet;

typedef nmatrix<double> TimePattern;

typedef DerivativeVector TestFunction;
typedef std::vector<TestFunction> FemFunction;
typedef std::map<std::string, FemFunction> FemData;

typedef nvector<double> CellFunction;
typedef std::map<std::string, CellFunction> CellData;

#define BASE_CLONEABLE(Type) virtual Type *createNew() const = 0;

#define CLONEABLE(Type)                                                        \
  virtual Type *createNew() const override { return new Type(*this); }
} // namespace Gascoigne

#endif
