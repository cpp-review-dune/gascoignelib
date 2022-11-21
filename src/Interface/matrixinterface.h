/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2009, 2011 by the Gascoigne 3D authors
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

#ifndef __matrixinterface_h
#define __matrixinterface_h

#include "compvector.h"
#include "entrymatrix.h"
#include "gascoigne.h"
#include "periodicdata.h"
#include "sparsestructureinterface.h"
#include "stencilinterface.h"
#include <string>

/*-------------------------------------------------------------*/

namespace Gascoigne {
class MatrixInterface
{
private:
protected:
public:
  MatrixInterface() {}
  virtual ~MatrixInterface() {}

  virtual std::string GetName() const = 0;

  virtual const StencilInterface* GetStencil() const = 0;
  virtual void ReInit(const SparseStructureInterface* S) = 0;
  virtual void AddMassWithDifferentStencil(
    [[maybe_unused]] const MatrixInterface* M,
    [[maybe_unused]] const TimePattern& TP,
    [[maybe_unused]] double s = 1.)
  {
    std::cerr << "\"MatrixInterface::AddMassWithDifferentStencil\" not written!"
              << std::endl;
    abort();
  }
  virtual void AddMassWithDifferentStencilJacobi(
    [[maybe_unused]] const MatrixInterface* M,
    [[maybe_unused]] const TimePattern& TP,
    [[maybe_unused]] double s = 1.)
  {
    std::cerr << "\"MatrixInterface::AddMassWithDifferentStencil\" not written!"
              << std::endl;
    abort();
  }
  virtual void copy_entries([[maybe_unused]] const MatrixInterface& S)
  {
    std::cerr << "\"MatrixInterface::copy_entries\" not written!" << std::endl;
    abort();
  }
  virtual void zero() = 0;
  virtual void transpose() = 0;

  virtual std::ostream& Write(std::ostream& os) const = 0;

  //
  /// for matrix assembling
  //
  typedef IntVector::const_iterator niiterator;

  virtual void entry([[maybe_unused]] nvector<int>::const_iterator start1,
                     [[maybe_unused]] nvector<int>::const_iterator stop1,
                     [[maybe_unused]] nvector<int>::const_iterator start2,
                     [[maybe_unused]] nvector<int>::const_iterator stop2,
                     [[maybe_unused]] const EntryMatrix& M,
                     [[maybe_unused]] double s = 1.)
  {
    std::cerr << "\"MatrixInterface::entry\" not written!" << std::endl;
    abort();
  }
  virtual void entry(niiterator start,
                     niiterator stop,
                     const EntryMatrix& M,
                     double s = 1.) = 0;
  virtual void entrydual([[maybe_unused]] niiterator start,
                         [[maybe_unused]] niiterator stop,
                         [[maybe_unused]] const EntryMatrix& M,
                         [[maybe_unused]] double s = 1.)
  {
    std::cerr << "\"MatrixInterface::entrydual\" not written!" << std::endl;
    abort();
  }

  //
  /// for hanging nodes
  //
  virtual void entry_diag(int i, const nmatrix<double>& M) = 0;

  //
  /// for boundary conditions
  //
  virtual void scale_diag([[maybe_unused]] int i,
                          [[maybe_unused]] const std::vector<int>& cv,
                          [[maybe_unused]] double s)
  {
    std::cerr << "\"MatrixInterface::scale_diag\" not written!" << std::endl;
    abort();
  }
  virtual void dirichlet([[maybe_unused]] int i,
                         [[maybe_unused]] const std::vector<int>& cv)
  {
    std::cerr << "\"MatrixInterface::dirichlet\" not written!" << std::endl;
    abort();
  }
  virtual void dirichlet_only_row(
    [[maybe_unused]] int i,
    [[maybe_unused]] const std::vector<int>& indices)
  {
    std::cerr << "\"MatrixInterface::dirichlet_only_row\" not written!"
              << std::endl;
    abort();
  }
  virtual void dirichlet_only_column(
    [[maybe_unused]] int i,
    [[maybe_unused]] const std::vector<int>& indices)
  {
    std::cerr << "\"MatrixInterface::dirichlet_only_column\" not written!"
              << std::endl;
    abort();
  }
  virtual void dirichlet_only_row_no_diag(
    [[maybe_unused]] int i,
    [[maybe_unused]] const std::vector<int>& indices)
  {
    std::cerr << "\"MatrixInterface::dirichlet_only_row_no_diag\" not written!"
              << std::endl;
    abort();
  }
  virtual void periodic(
    [[maybe_unused]] const std::map<int, int>& m_PeriodicPairs,
    [[maybe_unused]] const IntVector& iv_Components)
  {
    std::cerr << "\"MatrixInterface::periodic\" not written!" << std::endl;
    abort();
  }
  virtual void vmult([[maybe_unused]] GlobalVector& y,
                     [[maybe_unused]] const GlobalVector& x,
                     [[maybe_unused]] double s = 1.) const
  {
    std::cerr << "\"MatrixInterface::vmult\" not written!" << std::endl;
    abort();
  }
  virtual void vmult_transpose([[maybe_unused]] GlobalVector& y,
                               [[maybe_unused]] const GlobalVector& x,
                               [[maybe_unused]] double s = 1.) const
  {
    std::cerr << "\"MatrixInterface::vmult_tranpose\" not written!"
              << std::endl;
    abort();
  }
  virtual void vmult_time([[maybe_unused]] GlobalVector& y,
                          [[maybe_unused]] const GlobalVector& x,
                          [[maybe_unused]] const TimePattern& TP,
                          [[maybe_unused]] double s = 1.) const
  {
    std::cerr << "\"MatrixInterface::vmult_time\" not written!" << std::endl;
    abort();
  }

  /*-----------------------------------------------*/

  virtual void FillInterfaceList(
    [[maybe_unused]] const nvector<int>& elements,
    [[maybe_unused]] nvector<int>& start,
    [[maybe_unused]] nvector<MatrixEntryType>& values) const
  {
    std::cerr << "\"MatrixInterface::FillInterfaceList\" not written!"
              << std::endl;
    abort();
  }
  virtual void FurbishInterface(
    [[maybe_unused]] double d,
    [[maybe_unused]] const nvector<int>& elements,
    [[maybe_unused]] const nvector<int>& start,
    [[maybe_unused]] const nvector<MatrixEntryType>& values)
  {
    std::cerr << "\"MatrixInterface::FurbishInterface\" not written!"
              << std::endl;
    abort();
  }

  virtual void PrepareJacobi([[maybe_unused]] double s)
  {
    std::cerr << "\"MatrixInterface::PrepareJacobi\" not written!" << std::endl;
    abort();
  }
  virtual void Jacobi([[maybe_unused]] GlobalVector& x) const
  {
    std::cerr << "\"MatrixInterface::Jacobi\" not written!" << std::endl;
    abort();
  }

  /*-----------------------------------------------*/
};
} // namespace Gascoigne

/*-------------------------------------------------------------*/

#endif
