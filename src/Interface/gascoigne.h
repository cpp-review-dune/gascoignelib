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

#include "../Common/derivativevector.h"
#include "../Common/nmatrix.h"

namespace Gascoigne {

#ifdef __MATRIX_SINGLE_PRECISION__
typedef float MatrixEntryType;
#else
typedef double MatrixEntryType;
#endif

typedef int IndexType;      // data type for global indices
typedef int ShortIndexType; // data type for local indices

typedef nvector<double> GlobalParameterVector;
typedef nvector<double> LocalParameterVector;
typedef std::map<std::string, const GlobalParameterVector*> GlobalParameterData;
typedef std::map<std::string, LocalParameterVector> LocalParameterData;

typedef nvector<IndexType> IndexVector;
typedef std::set<IndexType> IndexSet;
typedef std::map<IndexType, IndexType> IndexMap;
typedef std::map<IndexType, IndexVector> VecMap;
typedef nvector<ShortIndexType> ShortIndexVector;
typedef nvector<int> IntVector;
typedef nvector<double> DoubleVector;
typedef nmatrix<double> DoubleMatrix;
typedef std::set<int> IntSet;

typedef nmatrix<double> TimePattern;

typedef DerivativeVector TestFunction;
typedef std::vector<TestFunction> FemFunction;
typedef std::map<std::string, FemFunction> FemData;

typedef nvector<double> CellFunction;
typedef std::map<std::string, CellFunction> CellData;

typedef std::array<IndexType, 2> EdgeVector;
typedef std::array<IndexType, 4> FaceVector;

#define BASE_CLONEABLE(Type) virtual Type* createNew() const = 0;

#define CLONEABLE(Type)                                                        \
  virtual Type* createNew() const override { return new Type(*this); }

/**
 * @brief Marker throwing error when not implemented
 *
 */
#define NOT_IMPLEMENTED                                                        \
  {                                                                            \
    throw std::runtime_error(std::string("Unimplemented funtion in ") +        \
                             std::string(__FILE__) + std::string(":") +        \
                             std::to_string(__LINE__) + std::string(" ") +     \
                             std::string(__FUNCTION__));                       \
  }
/**
 * @brief Marker throwing error when not implemented
 *
 */
#define ERROR(Error)                                                           \
  {                                                                            \
    throw std::runtime_error(std::string(Error) + std::string(" in ") +        \
                             std::string(__FILE__) + std::string(":") +        \
                             std::to_string(__LINE__) + std::string(" ") +     \
                             std::string(__FUNCTION__));                       \
  }

/**
 * @brief Marker to reminde of todos
 *
 */
#define TO_DO                                                                  \
  {                                                                            \
    printf("Todo in %s:%i %s\n", __FILE__, __LINE__, __FUNCTION__);            \
  }

/**
 * @brief Marker to reminde of todos
 *
 */
#define WARNING(Warning)                                                       \
  {                                                                            \
    printf(                                                                    \
      "WARNING: %s in %s:%i %s\n", Warning, __FILE__, __LINE__, __FUNCTION__); \
  }

/**
 * Helpfull marker for debuging purpose
 */
#define CHECK                                                                  \
  {                                                                            \
    std::cerr << "Past: " << __FILE__ << ":" << std::to_string(__LINE__)       \
              << " in " << __FUNCTION__ << std::endl;                          \
  }

} // namespace Gascoigne

#endif
