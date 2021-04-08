/*----------------------------   matrixagent.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __matrixagent_H
#define __matrixagent_H
/*----------------------------   matrixagent.h     ---------------------------*/

/**
 *
 * Copyright (C) 2004, 2005, 2020 by the Gascoigne 3D authors
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

#include "gascoigne.h"
#include "matrix.h"
#include "matrixinterface.h"

namespace Gascoigne {

/////////////////////////////////////////////
////
////@brief
////  stores the matrix data objects belonging to the Matrix labels
////
/////////////////////////////////////////////

class MatrixAgent : public std::map<Matrix, MatrixInterface *> {
public:
  typedef std::map<Matrix, MatrixInterface *>::const_iterator const_iterator;
  typedef std::map<Matrix, MatrixInterface *>::iterator iterator;

  //
  ////  Con(De)structor
  //
  MatrixAgent();
  ~MatrixAgent();

  void Register(const Matrix &mat);
  void Delete(Matrix &mat);

  MatrixInterface &operator()(const Matrix &g);
  const MatrixInterface &operator()(const Matrix &g) const;

  friend std::ostream &operator<<(std::ostream &os, const MatrixAgent &gva) {
    int i = 0, n = gva.size();
    os << "MatrixAgent: size=" << n << ", ";
    for (auto p = gva.begin(); p != gva.end(); p++, i++) {
      os << "Matrix(" << i << ")=('" << p->first << "'," << p->second << ")";
      if (i < n - 1)
        os << ", ";
      else
        os << ". ";
    }
    return os;
  }
};
} // namespace Gascoigne

/*----------------------------   matrixagent.h     ---------------------------*/
/* end of #ifndef __matrixagent_H */
#endif
/*----------------------------   matrixagent.h     ---------------------------*/
