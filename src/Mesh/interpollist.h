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


#ifndef __interpollist_h
#define __interpollist_h

#include  <vector.h>
#include  "fixarray.h" 

/*---------------------------------------------------*/

namespace Gascoigne
{
template <int N>
class InterpolElement : public std::array<int,N>
{
 public:

  int nv;

  InterpolElement() : std::array<int,N>(), nv(-1) {}
  InterpolElement(const InterpolElement& i) : std::array<int,N>(i), 
    nv(i.nv) {}

  InterpolElement(const std::array<int,N>& f, int n) : 
    std::array<int,N>(f), nv(n) {}
};

/*---------------------------------------------------*/

template <int N>
class InterpolationList : public std::vector<InterpolElement<N> >
{
  public :
    
  int  newvertex(int i)        const { return (*this)[i].nv; }
  int  oldvertex(int i, int j) const { return (*this)[i][j]; }

  void newentry(int nv, const std::array<int,N>& w)
    {
      InterpolElement<N> I(w,nv);
      push_back(I);
    }
};
}

#endif
