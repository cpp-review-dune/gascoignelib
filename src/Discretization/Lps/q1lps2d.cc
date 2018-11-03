/**
*
* Copyright (C) 2004, 2009 by the Gascoigne 3D authors
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


#include  "q1lps2d.h"

namespace Gascoigne
{
/* ----------------------------------------- */

Q1Lps2d:: ~Q1Lps2d()
{
  if (S) delete S;
  S=0;
}

/* ----------------------------------------- */  


void Q1Lps2d::BasicInit(const ParamFile* paramfile)
{
  Q12d::BasicInit(paramfile);
  S = new Q1LpsStab2d;
  S ->BasicInit(paramfile,HN);
}

/* ----------------------------------------- */

void Q1Lps2d::ReInit(const GascoigneMesh* M)
{
  Q12d::ReInit(M);
  S   ->ReInit(M);
}

/* ----------------------------------------- */

void Q1Lps2d::Structure(SparseStructureInterface* SI) const
{
  S->Structure(SI);
}


/* ----------------------------------------- */

void Q1Lps2d::StabForm(GlobalVector& f, const GlobalVector& u, const ProblemDescriptorInterface& PD, double d) const
{
  assert(0);
  
  //  S->Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps2d::Form(GlobalVector& f, const GlobalVector& u, const ProblemDescriptorInterface& PD, double d) const
{
  Q12d::Form(f,u,PD,d);
  S   ->Form(f,u,PD,d);
}

/* ----------------------------------------- */

void Q1Lps2d::Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface& PD, double d) const
{
  Q12d::Matrix(A,u,PD,d);
  S   ->Matrix(A,u,PD,d);
}
}
