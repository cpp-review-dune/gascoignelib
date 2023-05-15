/*----------------------------   navierstokesproblem.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __navierstokesproblem_H
#define __navierstokesproblem_H
/*----------------------------   navierstokesproblem.h     ---------------------------*/


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


#include <navierstokes.h>     // see Gascoigne/src/Problem/navierstokes.h
#include "navierstokesimex.h"
#include "dirichletdata.h"
#include "paramfile.h"
#include "problemdescriptorbase.h"


/*-----------------------------------------*/

namespace Gascoigne
{
  
  /**
   * \brief Defines the DirichletData of the problem.
   *
   **/
  class BenchDD : public DirichletData
  {
  protected:

    double vmax;   // maximum velocity at inflow
    
  public:
    /** 
     * \brief Constructor
     *
     * We must always pass the parameter file to the constructor of
     * the base class. Otherwise, the correct boundary information
     * is not taken from the param file
     */
    BenchDD(const ParamFile& pf) : DirichletData(pf)
    {
      DataFormatHandler DFH;
      DFH.insert("vmax", &vmax, 0.);
      FileScanner FS(DFH, pf, "Equation");
    }


    /// @brief Dirichlet Data in 2d
    void operator()(DoubleVector &b, const Vertex2d &v, int color) const
    {
      b.zero();
      if (color==8) // 8 is the color of the inflow boundary in the mesh nsbench
	b[1] = vmax * v.y()*(0.41-v.y()) / 0.205/0.205;
    }
    /// @brief Dirichlet Data in 3d
    void operator()(DoubleVector &b, const Vertex3d &v, int color) const
    {
      b.zero();
      if (color==8) // 8 is the color of the inflow boundary in the mesh nsbench
	b[1] = vmax 
	  * v.y()*(0.41-v.y()) / 0.205/0.205
	  * v.z()*(0.41-v.z()) / 0.205/0.205;
    }
  };



  /**
   * \brief Description of the 2d-1 Benchmark problem
   *
   */
  template<int DIM>
  class NavierStokesStationary : public ProblemDescriptorBase
  {
    NavierStokesData data;
    
  public:

    const NavierStokesData& GetData() const  { return data; }
    
    // Setting up the problem and the data
    void BasicInit(const ParamFile &pf)
    {
      data.BasicInit(pf);
      
      GetParamFile()     = pf;
      
      GetEquationPointer()      = new Gascoigne::NavierStokesLps<DIM>(data);
      GetDirichletDataPointer() = new BenchDD(pf);
      ProblemDescriptorBase::BasicInit(pf);
    }
  };

  /**
   * \brief Description of the instationary NS Benchmark problem
   *
   */
  template<int DIM>
  class NavierStokesInstationary : public ProblemDescriptorBase
  {
    NavierStokesData data;
    
  public:

    const NavierStokesData& GetData() const  { return data; }
    
    // Setting up the problem and the data
    void BasicInit(const ParamFile &pf)
    {
      data.BasicInit(pf);
      
      GetParamFile()            = pf;
      GetEquationPointer()      = new Gascoigne::NavierStokesLpsTime<DIM>(data);
      GetDirichletDataPointer() = new BenchDD(pf);
      ProblemDescriptorBase::BasicInit(pf);
    }
  };

  /**
   * \brief Description of the instationary NS Benchmark problem using IMEX schemes (BDF2)
   *
   */
  template<int DIM>
  class NavierStokesInstationaryIMEX : public ProblemDescriptorBase
  {
    NavierStokesData data;
    
  public:

    const NavierStokesData& GetData() const  { return data; }
    
    // Setting up the problem and the data
    void BasicInit(const ParamFile &pf)
    {
      data.BasicInit(pf);
      
      GetParamFile()            = pf;
      GetEquationPointer()      = new Gascoigne::NavierStokesLpsTimeIMEX<DIM>(data);
      GetDirichletDataPointer() = new BenchDD(pf);
      ProblemDescriptorBase::BasicInit(pf);
    }
  };




}

/*----------------------------   navierstokesproblem.h     ---------------------------*/
/* end of #ifndef __navierstokesproblem_H */
#endif
/*----------------------------   navierstokesproblem.h     ---------------------------*/
