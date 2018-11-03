/*----------------------------   navierstokes.h ---------------------------*/
/*      $Id:$                 */
#ifndef __navierstokes_H
#define __navierstokes_H
/*----------------------------   navierstokes.h ---------------------------*/


/**
 *
 * Copyright (C) 2004,2018 by the Gascoigne 3D authors
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


#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"


/*-----------------------------------------*/

namespace Gascoigne
{

  class NavierStokesData
  {
  public:
    double visc, alpha0, dt;

    void BasicInit(const ParamFile *pf);
  };

  template <int DIM>
  class NavierStokes : public virtual Equation
  {

  protected:
    NavierStokesData data;

  public:
    NavierStokes<DIM>(const NavierStokesData &PD)
        : data(PD)
    {
    }

    std::string GetName() const
    { return "NavierStokes"; }
    int GetNcomp() const
    {
      return DIM + 1;
    }

    //
    // Semilinear Form
    //

    void
    Form(VectorIterator b, const FemFunction &U, const TestFunction &N) const;
    void Matrix(EntryMatrix &A,
                const FemFunction &U,
                const TestFunction &M,
                const TestFunction &N) const;
  };


  template <int DIM>
  class NavierStokesLps : public virtual NavierStokes<DIM>, public LpsEquation
  {
    mutable double _h, _alpha;

  public:
    NavierStokesLps<DIM>(const NavierStokesData &PD)
        : NavierStokes<DIM>(PD)
    {
    }
    std::string GetName() const
    { return "NavierStokesLps"; }
    

    void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const;
    void StabForm(VectorIterator b,
                  const FemFunction &U,
                  const FemFunction &UP,
                  const TestFunction &NP) const;
    void StabMatrix(EntryMatrix &A,
                    const FemFunction &U,
                    const TestFunction &Np,
                    const TestFunction &Mp) const;
  };

  template <int DIM>
  class NavierStokesLpsTime : public virtual NavierStokesLps<DIM>
  {
    mutable FemFunction *OLD;
    

  public:
    NavierStokesLpsTime<DIM>(const NavierStokesData &PD)
      : NavierStokes<DIM>(PD), NavierStokesLps<DIM>(PD)
    {
    }
    std::string GetName() const
    { return "NavierStokesTimeLps"; }
    

    void SetFemData(FemData& q) const
    {
      assert(q.find("OLD")!=q.end());
      OLD = &q["OLD"];
    }

    void
    Form(VectorIterator b, const FemFunction &U, const TestFunction &N) const
    {      
      NavierStokesLps<DIM>::Form(b, U, N);
      
      for (int i = 0; i < DIM; ++i)
        b[i+1] += (U[i + 1].m() - (*OLD)[i+1].m()) / NavierStokes<DIM>::data.dt * N.m();
    }

    void Matrix(EntryMatrix &A,
                const FemFunction &U,
                const TestFunction &M,
                const TestFunction &N) const
    {
      NavierStokesLps<DIM>::Matrix(A, U, M, N);

      for (int i = 0; i < DIM; ++i)
        A(i+1, i+1) += M.m() / NavierStokes<DIM>::data.dt * N.m();
    }
  };

#define NavierStokes2d NavierStokes<2>
#define NavierStokes3d NavierStokes<3>
#define NavierStokesLps2d NavierStokesLps<2>
#define NavierStokesLps3d NavierStokesLps<3>
  
#define NavierStokesLpsTime2d NavierStokesLpsTime<2>
#define NavierStokesLpsTime3d NavierStokesLpsTime<3>

} // namespace Gascoigne


/*----------------------------   navierstokes.h ---------------------------*/
/* end of #ifndef __navierstokes_H */
#endif
/*----------------------------   navierstokes.h ---------------------------*/
