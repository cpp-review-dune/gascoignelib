/**
*
* Copyright (C) 2004, 2008, 2011 by the Gascoigne 3D authors
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


#include "gascoignehash.h"
#include "hanglist.h"
#include <cassert>
#include "stlio.h"

using namespace std;

#define HANGMAP hash_map<EdgeArray<N>, Hang, EdgeHash<N>>


/*********************************************************************/

namespace Gascoigne
{
  template <int N>
  ostream &operator<<(ostream &s, const HangList<N> &A)
  {
    s << A.size() << " hangs" << endl;
    auto p = A.begin();

    for (p = A.begin(); p != A.end(); p++)
    {
      s << p->first << "-> " << p->second;
    }
    return s;
  }

  /*********************************************************************/

  template <int N>
  istream &operator>>(istream &s, HangList<N> &A)
  {
    int n;
    string symbol;

    s >> n >> symbol;

    assert(symbol == "hangs");

    std::array<int, N> ev;
    Hang info;
    for (int i = 0; i < n; i++)
    {
      s >> ev >> symbol >> info;
      A.insert(make_pair(ev, info));
    }

    return s;
  }

  /*********************************************************************/

  template ostream &operator<<(ostream &s, const HangList<2> &A);
  template ostream &operator<<(ostream &s, const HangList<4> &A);
  template istream &operator>>(istream &s, HangList<2> &A);
  template istream &operator>>(istream &s, HangList<4> &A);
}

#undef HANGMAP
