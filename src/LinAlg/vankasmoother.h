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

#include "dofhandler.h"
#include "iluinterface.h"
#include "sparseblockmatrix.h"
#include "fmatrixblock.h"
#include "matrixentrytype.h"

// we use Eigen to store and invert the local matrices
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wold-style-cast"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/OrderingMethods>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#pragma clang diagnostic pop


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
namespace Gascoigne
{
  class VankaSmoother : public virtual IluInterface
  {
    mutable const DofHandlerBase *_dofhandler;
    int _ncomp, _sizeofpatch;

  protected:
    typedef Eigen::Matrix<MatrixEntryType, Eigen::Dynamic, Eigen::Dynamic> VankaMatrix;
    typedef Eigen::Matrix<MatrixEntryType, Eigen::Dynamic, 1>              VankaVector;

    std::vector<std::vector<int>> _patchlist;
    mutable std::vector<Eigen::PartialPivLU<VankaMatrix> > _lu;
    std::vector<std::vector<int> >  Coloring;
  public:
    //////////////////// Constructor & Co
    VankaSmoother() : _dofhandler(NULL), _ncomp(-1), _sizeofpatch(-1)
    {
    }
    ~VankaSmoother() {}

    void SetDofHandler(const DofHandlerBase* dh) const
    {
      _dofhandler = dh;
    }

    std::string GetName() const
    {
      return "VankaSmoother";
    }

    //////////////////// Access
    int n() const
    {
      assert(0);
    }
    void ReInit(const SparseStructureInterface *A)
    {
      // nothing to be done, Vanka smoother does not depend on the stencil
    }

    //////////////////// Construction
    void ConstructStructure(const IntVector &perm, const MatrixInterface &A);
    void zero()
    {
      // not necessary. entries will be copied
    }
    void modify(int c, double s)
    {
      // not necessary
    }
    void compute_ilu ()
    {
      // directly done in ConstructStructure
    }

    //////////////////// Solve
    void solve(GlobalVector& x) const;


    template<int NCOMP>
    void copy_entries_sparseblockmatrix(const SparseBlockMatrix<FMatrixBlock<NCOMP> >& A);
    void copy_entries(const MatrixInterface& A);


    void ElementColoring(int degree)
    {
      // convert mesh to graph
      std::vector<std::vector<int> > node2patch(_dofhandler->nnodes());

      for (int p=0;p<_patchlist.size();++p)
	{
	  const auto &iop = _patchlist[p];
	  for (auto it : iop)
	    node2patch[it].push_back(p);
	}
      std::vector<std::vector<int> > adj(_patchlist.size());
      for (int p=0;p<_patchlist.size();++p)
	{
	  std::set<int> neigh;
	  const auto &iop =_patchlist[p];
	  for (auto it : iop)
	    neigh.insert(node2patch[it].begin(),node2patch[it].end());
	  for (auto it : neigh)
	    if (it!=p)
	      adj[p].push_back(it);
	}
      /////////
      std::vector<int> patch2color(_patchlist.size(),0);
      bool done = true;
      int color = 0;
      //faerben aller patches
      do
	{
	  done = true;
	  for (int p=0;p<_patchlist.size();++p)
	    if ((patch2color[p]==color) ) // not given a color and fluid
	      {
		done = false;
		// block neighbors
		for (auto n : adj[p])
		  {
		    if ((patch2color[n] == color) ) // neighbor not set, block
		      patch2color[n] = color+1;
		  }
	      }
	  ++color;
	} while (!done);

      Coloring.clear();
      Coloring.resize(color-1);
      for (int p=0;p<_patchlist.size();++p)
	{
	  int col = patch2color[p];
	  assert(col<Coloring.size());
	  Coloring[col].push_back(p);
	}
    }
    int NumberofColors() const
    {
      return Coloring.size();
    }

    const std::vector<int>&   elementswithcolor(int col) const
    {
      return Coloring[col];
    }
  };
} // namespace Gascoigne


/*----------------------------   vank_smoother.h ---------------------------*/
/* end of #ifndef __vanksmoother_H */
#endif
/*----------------------------   vank_smoother.h ---------------------------*/
