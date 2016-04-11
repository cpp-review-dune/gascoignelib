/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "myeq.h"
#include  "filescanner.h"


extern double __DT;


/*-----------------------------------------*/

namespace Gascoigne
{


MyEQ::~MyEQ()
{
}

/*-----------------------------------------*/

MyEQ::MyEQ() : Equation()
{
  _visc = 0.01;
}

/*-----------------------------------------*/

MyEQ::MyEQ(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 1.);
  DFH.insert("lps",&_lps0, 0.);

  FileScanner FS(DFH, pf, "Equation");
  assert(_lps0>0);
}

/*-----------------------------------------*/

void MyEQ::OperatorStrong(DoubleVector& b, const FemFunction& U) const
{
}

/*-----------------------------------------*/

double MyEQ::Laplace(const TestFunction& U, const TestFunction& N) const
{
  return U.x()*N.x() + U.y()*N.y();
}

/*-----------------------------------------*/

double MyEQ::Convection(const FemFunction& U, const TestFunction& N) const
{
  return U[1].m()*N.x() + U[2].m()*N.y();
}

/*-----------------------------------------*/

double MyEQ::Divergence(const FemFunction& U) const
{
  return U[1].x() + U[2].y();
}
 
/*-----------------------------------------*/

void MyEQ::SetTimePattern(TimePattern& P) const
{
}

/*-----------------------------------------*/

void MyEQ::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  ////////////// Continuity ////////////////////////////////////////////////

  b[0] += Divergence(U) * N.m();

  ////////////// Momentum ////////////////////////////////////////////////

  b[1] += (U[1].m()-(*old)[1].m())*N.m()/__DT;
  b[2] += (U[2].m()-(*old)[2].m())*N.m()/__DT;
  

  b[1] -= U[0].m()*N.x();
  b[2] -= U[0].m()*N.y();
	  
  b[1] += Convection(U,U[1]) * N.m();
  b[2] += Convection(U,U[2]) * N.m();

  // viscous terms
  b[1] += _visc * Laplace(U[1],N);
  b[2] += _visc * Laplace(U[2],N);
}

/*-----------------------------------------*/

void MyEQ::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  double MN = M.m()*N.m();
  double Mx = M.x()*N.m();
  double My = M.y()*N.m();
  double laplace = Laplace(M,N);

  A(1,1) += MN/__DT;
  A(2,2) += MN/__DT;
  
  
  ////////////// Continuity ////////////////////////////////////////////////

  A(0,1) += Mx;
  A(0,2) += My;
  
  ////////////// Momentum ////////////////////////////////////////////////

  A(1,0) -= M.m()*N.x();
  A(2,0) -= M.m()*N.y();
  
  double cl = Convection(U,M) * N.m() + _visc*laplace;
  
  A(1,1) += cl;
  A(2,2) += cl;
  
  A(1,1) += U[1].x()*MN;
  A(2,2) += U[2].y()*MN;
  A(1,2) += U[1].y()*MN;
  A(2,1) += U[2].x()*MN;
}
}
/*-----------------------------------------*/
