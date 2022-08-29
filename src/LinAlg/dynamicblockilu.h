/**
 *
 * Copyright (C) 2008, 2011 by the Gascoigne 3D authors
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

#ifndef __dynamicblockilu_h
#define __dynamicblockilu_h

#include "dynamicblockmatrix.h"
#include "iluinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne {
template<class B>
class DynamicBlockIlu
  : public virtual IluInterface
  , public DynamicBlockMatrix<B>
{
protected:
  typedef typename std::list<IndexType>::const_iterator const_citerator;
  typedef typename std::list<IndexType>::iterator citerator;
  typedef typename std::list<B>::const_iterator const_viterator;
  typedef typename std::list<B>::iterator viterator;

  nvector<IndexType> p, q;
  GlobalVector* yp;

  void backward() const;
  void forward() const;
  virtual void hin(const GlobalVector& x) const;
  virtual void her(GlobalVector& x) const;

  IndexType n() const { return DynamicBlockMatrix<B>::n(); };

public:
  DynamicBlockIlu();
  DynamicBlockIlu(const DynamicBlockIlu<B>& I);
  ~DynamicBlockIlu();

  viterator vdiag(IndexType i) { return DynamicBlockMatrix<B>::vdiag(i); }
  const_viterator vdiag(IndexType i) const
  {
    return DynamicBlockMatrix<B>::vdiag(i);
  }
  citerator cdiag(IndexType i) { return DynamicBlockMatrix<B>::cdiag(i); }
  const_citerator cdiag(IndexType i) const
  {
    return DynamicBlockMatrix<B>::cdiag(i);
  }
  const_citerator cstart(IndexType i) const
  {
    return DynamicBlockMatrix<B>::cstart(i);
  }
  const_citerator cstop(IndexType i) const
  {
    return DynamicBlockMatrix<B>::cstop(i);
  }
  citerator cstart(IndexType i) { return DynamicBlockMatrix<B>::cstart(i); }
  citerator cstop(IndexType i) { return DynamicBlockMatrix<B>::cstop(i); }
  const_viterator vstart(IndexType i) const
  {
    return DynamicBlockMatrix<B>::vstart(i);
  }
  const_viterator vstop(IndexType i) const
  {
    return DynamicBlockMatrix<B>::vstop(i);
  }
  viterator vstart(IndexType i) { return DynamicBlockMatrix<B>::vstart(i); }
  viterator vstop(IndexType i) { return DynamicBlockMatrix<B>::vstop(i); }

  std::string GetName() const { return "DynamicBlockIlu"; }

  nvector<IndexType>& GetP() { return p; }
  nvector<IndexType>& GetQ() { return q; }
  const nvector<IndexType>& GetP() const { return p; }
  const nvector<IndexType>& GetQ() const { return q; }

  void modify(IndexType c, double s);
  void zero() { DynamicBlockMatrix<B>::zero(); }

  void compute_ilu();
  void ReInit(const SparseStructureInterface* SI);
  void ConstructStructure(const nvector<IndexType>& perm,
                          const MatrixInterface& A);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  void copy_entries(const MatrixInterface* A);
#pragma GCC diagnostic pop
  void solve(GlobalVector& x) const;
  void solvetrans(GlobalVector& x) const
  {
    std::cerr << "\"DynamicBlockIlu::solvetrans\" not written!" << std::endl;
    abort();
  };
};
} // namespace Gascoigne

#endif
