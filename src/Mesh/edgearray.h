/**
 *
 * Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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

#ifndef __EdgeArray_h
#define __EdgeArray_h

#include <array>
#include <iostream>

#include "../Interface/gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne {
template<size_t N>
class EdgeArray : public std::array<IndexType, N>
{
public:
  EdgeArray(const std::array<IndexType, N>& e);

  bool operator==(const std::array<IndexType, N>& e) const;

  IndexType sum() const;

  void BinWrite(std::ostream& s) const
  {
    int sizeT = sizeof(int);
    s.write(reinterpret_cast<const char*>(this->data()), sizeT * N);
  }

  void BinRead(std::istream& s)
  {
    int sizeT = sizeof(int);
    s.read(reinterpret_cast<char*>(this->data()), sizeT * N);
  }
};
} // namespace Gascoigne

/*------------------------------------------------------*/

#endif
