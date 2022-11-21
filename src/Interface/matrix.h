/*----------------------------   matrix.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __matrix_H
#define __matrix_H
/*----------------------------   matrix.h     ---------------------------*/

/**
 *
 * Copyright (C) 2020 by the Gascoigne 3D authors
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

#include <cassert>
#include <iostream>
#include <string>

namespace Gascoigne {
/////////////////////////////////////////////
////
////@brief
////  Matrix is just a placeholder for any kind of (large) matrices that act on
////  all dof's of a discretization. This class takes the same role as the
////  Vector (to be renamed as Vector).
////  The Matrix class does not have any real data. The connection to the data
////  is by means of the MatrixAgent class.
////  However, the Matrix class stores some extra information such as
////
////  -
/////////////////////////////////////////////

class Matrix : public std::string
{
private:
protected:
public:
  Matrix()
  {
    std::cerr << "Matrix: no Constructor without indicating the label!"
              << std::endl;
    abort();
  }

  Matrix(const std::string& label)
    : std::string(label)
  {
  }
  Matrix(const Matrix& mat)
    : std::string(mat)
  {
  }
  virtual ~Matrix() {}
};
} // namespace Gascoigne

/*----------------------------   matrix.h     ---------------------------*/
/* end of #ifndef __matrix_H */
#endif
/*----------------------------   matrix.h     ---------------------------*/
