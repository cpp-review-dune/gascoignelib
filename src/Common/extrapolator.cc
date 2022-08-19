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

#include "extrapolator.h"
#include <iomanip>

using namespace std;

namespace Gascoigne {

/*-----------------------------------------*/

Extrapolator::Extrapolator() {}
Extrapolator::~Extrapolator() {}

/*-----------------------------------------*/

void
Extrapolator::Print()
{
  size_t n = order.size();
  cout << "extra\t";

  cout.precision(12);
  for (int i = 0; i < n; i++) {
    cout << fixed << std::setw(16) << valextra[i];
  }
  cout << endl << "[rate]\t";
  for (int i = 0; i < n; i++) {
    cout.precision(4);
    cout << setw(8) << "[";
    cout << setw(4) << order[i] << "]";
    cout.precision(12);
  }
  cout << endl;
}

/*-----------------------------------------*/

void
Extrapolator::NewValues(const DoubleVector& J)
{
  size_t n = J.size();
  if (!vals.size()) {
    valextra.resize(n);
    order.resize(n);
    valextra = 0.;
    order = -1.;
  }
  vals.push_back(J);
  size_t nv = vals.size();
  if (nv >= 3) {
    for (int i = 0; i < n; i++) {
      double j0 = vals[nv - 3][i];
      double j1 = vals[nv - 2][i];
      double j2 = vals[nv - 1][i];

      double d0 = j0 - j1;
      double d1 = j1 - j2;
      double alpha = d0 / d1;
      double c = alpha * d0 / (alpha - 1.);
      double j = j0 - c;

      order[i] = alpha;
      valextra[i] = j;
    }
  }
}
} // namespace Gascoigne
