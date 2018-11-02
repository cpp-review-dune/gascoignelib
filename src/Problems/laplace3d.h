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


#ifndef  __Laplace3d_h
#define  __Laplace3d_h

#include  "laplace2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{

  class LaplaceData3d : virtual public LaplaceData2d
  {
  public:
    double betax, betay, betaz;

    LaplaceData3d() : LaplaceData2d()
    {
    }
    
    LaplaceData3d(const ParamFile* pf) : LaplaceData2d(pf)
    {
      DataFormatHandler DFH;
      DFH.insert("betax",&betax,0.);
      DFH.insert("betay",&betay,0.);
      DFH.insert("betaz",&betaz,0.);
      FileScanner FS(DFH, pf, "Equation");
    }
    
  };
  
class Laplace3d : public Laplace2d
{
  double betax, betay, betaz;

public:

  Laplace3d(const LaplaceData3d& data) : Laplace2d(data)
  {
    betax = data.betax;
    betay = data.betay;
    betaz = data.betaz;
  }
  

  std::string GetName() const { return "Laplace3d";}

  //
  // Semilinear Form
  //

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& D, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;

  void OperatorStrong(DoubleVector& b, const FemFunction& U)const;
};
}

#endif
