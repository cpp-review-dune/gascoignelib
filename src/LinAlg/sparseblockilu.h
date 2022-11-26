/**
 *
 * Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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

#ifndef __sparseblockilu_h
#define __sparseblockilu_h

#include "../Interface/iluinterface.h"

#include "sparseblockmatrix.h"

/*-------------------------------------------------------------*/

namespace Gascoigne {
template<class B>
class SparseBlockIlu
  : public virtual IluInterface
  , public SparseBlockMatrix<B>
{
protected:
  nvector<IndexType> p, q;
  GlobalVector* yp;

  void backward() const;
  void forward() const;
  virtual void hin(const GlobalVector& x) const;
  virtual void her(GlobalVector& x) const;

  IndexType n() const { return SparseBlockMatrix<B>::US.n(); };
  const IndexType& start(IndexType i) const
  {
    return SparseBlockMatrix<B>::US.start(i);
  };
  const IndexType& stop(IndexType i) const
  {
    return SparseBlockMatrix<B>::US.stop(i);
  };
  const IndexType& col(IndexType pos) const
  {
    return SparseBlockMatrix<B>::US.col(pos);
  };
  const IndexType& diag(IndexType i) const
  {
    return SparseBlockMatrix<B>::US.diag(i);
  };

public:
  SparseBlockIlu();
  SparseBlockIlu(const SparseBlockIlu<B>& I);
  ~SparseBlockIlu();

  std::string GetName() const { return "SparseBlockIlu"; }

  void dirichletILU(IndexType i, const std::vector<IndexType>& cv);

  nvector<IndexType>& GetP() { return p; }
  nvector<IndexType>& GetQ() { return q; }
  const nvector<IndexType>& GetP() const { return p; }
  const nvector<IndexType>& GetQ() const { return q; }

  void modify(IndexType c, double s);
  void zero() { SparseBlockMatrix<B>::zero(); }

  void compute_ilu();
  void ReInit(const SparseStructureInterface* SI);
  void ConstructStructure(const nvector<IndexType>& perm,
                          const MatrixInterface& A);

  void copy_entries(const MatrixInterface& A);

  void solve(GlobalVector& x) const;
  void solvetrans(GlobalVector& x) const
  {
    std::cerr << "\"SparseBlockIlu::solvetrans\" not written!" << std::endl;
    abort();
  }
  std::ostream& Write(std::ostream& s) const;
};
} // namespace Gascoigne

#endif
