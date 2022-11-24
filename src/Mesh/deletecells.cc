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

#include "deletecells.h"
#include "boundarycell.h"
#include "cell.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {
template<class C>
void
delete_cells(const IndexSet& coarselist,
             vector<C>& liste,
             const IndexVector& co2n,
             const IndexVector& vo2n)
{
  for (unsigned oi = 0; oi < co2n.size(); ++oi) {
    int ni = co2n[oi];
    if (ni >= 0) {
      C q(liste[oi]);

      for (unsigned i = 0; i < q.vertex().size(); ++i) {
        q.vertex(i) = vo2n[q.vertex(i)];
        if (q.vertex(i) == -1) {
          cerr << "Vertex invalid in " << oi << " " << ni << endl;
          // cerr << vo2n[liste[oi].vertex(i)]<<endl;
          // cerr << q.vertex();
          abort();
        }
      }
      IndexSet::iterator p = coarselist.find(oi);
      if (p != coarselist.end()) {
        q.childs().resize(0);
      }
      int qf = -1;
      if (q.father() != -1)
        qf = co2n[q.father()];
      q.father() = qf;
      if (q.sleep()) {
        for (int i = 0; i < q.nchilds(); ++i) {
          q.child(i) = co2n[q.child(i)];
        }
      }
      liste[ni] = q;
    }
  }
}

/*---------------------------------------------------*/

template void
delete_cells<Quad>(const IndexSet&,
                   vector<Quad>&,
                   const IndexVector&,
                   const IndexVector&);
template void
delete_cells<Hex>(const IndexSet&,
                  vector<Hex>&,
                  const IndexVector&,
                  const IndexVector&);

template void
delete_cells<BoundaryCell<2>>(const IndexSet&,
                              vector<BoundaryCell<2>>&,
                              const IndexVector&,
                              const IndexVector&);

template void
delete_cells<BoundaryCell<4>>(const IndexSet&,
                              vector<BoundaryCell<4>>&,
                              const IndexVector&,
                              const IndexVector&);
} // namespace Gascoigne
