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

#ifndef __visudata_h
#define __visudata_h

#include "vertex.h"

/***************************************************************/

namespace Gascoigne {
class VisuData
{
private:
protected:
public:
  VisuData() {}
  virtual ~VisuData() {}

  virtual int visucomp() const { return 0; }
  virtual int visun() const { return 0; }
  virtual double visudata(int i, int c) const
  {
    std::cerr << "\"VisuData::visudata\" not written!" << std::endl;
    abort();
  }
  virtual double visudata2(int i, int c, const Vertex2d& v) const
  {
    return visudata(i, c);
  }
  virtual double visudata2(int i, int c, const Vertex3d& v) const
  {
    return visudata(i, c);
  }
};
} // namespace Gascoigne

/***************************************************************/

#endif
