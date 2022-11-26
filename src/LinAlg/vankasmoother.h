/*----------------------------   vank_smoother.h ---------------------------*/
/*      $Id:$                 */
#ifndef __vanksmoother_H
#define __vanksmoother_H
/*----------------------------   vank_smoother.h ---------------------------*/

/**
 *
 * Copyright (C) 2018 by the Gascoigne 3D authors
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

#include "../DofHandler/dofhandler.h"
#include "../Interface/gascoigne.h"
#include "../Interface/iluinterface.h"

#include "fmatrixblock.h"
#include "sparseblockmatrix.h"

// we use Eigen to store and invert the local matrices
//#pragma clang diagnostic push
//#pragma clang diagnostic ignored "-Wunused-variable"
//#pragma clang diagnostic ignored "-Wold-style-cast"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/OrderingMethods>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
//#pragma clang diagnostic pop
#pragma GCC diagnostic pop

/*-------------------------------------------------------------*/

/**
 * Vanka smoother:
 *
 * the dofs are clusted in patches. The variable _patchlist
 * collects indices that belong to a certain patch.
 * - The restriction of the matrix to these patches is locally copied
 *   its inverse is stored as LU-decomposition in _lu
 *
 * Solve performs an Jacobi-iteration with averaging on the common
 * nodes for the patches.
 **/
namespace Gascoigne {
class VankaSmoother : public virtual IluInterface
{

protected:
  typedef Eigen::Matrix<MatrixEntryType, Eigen::Dynamic, Eigen::Dynamic>
    VankaMatrix;
  typedef Eigen::Matrix<MatrixEntryType, Eigen::Dynamic, 1> VankaVector;

  mutable const DofHandlerBase* _dofhandler;
  IndexType _ncomp, _sizeofpatch;

  std::vector<std::vector<IndexType>> _patchlist;
  mutable std::vector<Eigen::PartialPivLU<VankaMatrix>> _lu;

public:
  //////////////////// Constructor & Co
  VankaSmoother()
    : _dofhandler(NULL)
    , _ncomp(-1)
    , _sizeofpatch(-1)
  {}
  ~VankaSmoother() {}

  void SetDofHandler(const DofHandlerBase* dh) const { _dofhandler = dh; }

  std::string GetName() const { return "VankaSmoother"; }

  //////////////////// Access
  IndexType n() const
  {
    assert(0);
    return 0;
  }
  void ReInit(const SparseStructureInterface* A)
  {
    // nothing to be done, Vanka smoother does not depend on the stencil
  }

  //////////////////// Construction
  void ConstructStructure(const IndexVector& perm, const MatrixInterface& A);
  void zero()
  {
    // not necessary. entries will be copied
  }
  void modify(IndexType c, double s)
  {
    // not necessary
  }
  void compute_ilu()
  {
    // directly done in ConstructStructure
  }

  //////////////////// Solve
  void solve(GlobalVector& x) const;

  template<int NCOMP>
  void copy_entries_sparseblockmatrix(
    const SparseBlockMatrix<FMatrixBlock<NCOMP>>& A);
  void copy_entries(const MatrixInterface& A);
};
} // namespace Gascoigne

/*----------------------------   vank_smoother.h ---------------------------*/
/* end of #ifndef __vanksmoother_H */
#endif
/*----------------------------   vank_smoother.h ---------------------------*/
