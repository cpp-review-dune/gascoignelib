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

#ifndef __dirichletdata_h
#define __dirichletdata_h

#include "application.h"
#include "filescanner.h"
#include "nvector.h"
#include "vertex.h"
#include <set>
#include <string>

/*-----------------------------------------*/

namespace Gascoigne {

//////////////////////////////////////////////
///
///@brief
/// Interface class for Dirichlet Boundary Conditions

/// void operator()(Vector& b, const Vertex2d& v, int col)
/// gets the coordinate v and color of boundarypart "col" and
/// sets the values of b. b is a vector of length ncomp
///
//////////////////////////////////////////////

class DirichletData : public virtual Application
{
private:
protected:
  std::set<int> colors; // colors where Dirichlet data is given
  std::map<int, IntVector> comp_on_color; // components for each color

public:
  DirichletData(const ParamFile& pf)
  {
    DataFormatHandler DF;
    DF.insert("dirichlet", &colors);
    DF.insert("dirichletcomp", &comp_on_color);
    FileScanner FS(DF, pf, "BoundaryManager");
  }

  DirichletData()
  {
    std::cerr << "Warning: DirichletData without colors and comp_on_color"
              << std::endl;
  }
  virtual ~DirichletData() {}

  virtual void operator()(DoubleVector& b, const Vertex2d& v, int col) const
  {
    std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
    abort();
  }

  virtual void operator()(DoubleVector& b, const Vertex3d& v, int col) const
  {
    std::cerr << "\"DirichletData::operator()\" not written!" << std::endl;
    abort();
  }

  virtual std::set<int> preferred_colors() const { return std::set<int>(); }

  virtual const std::set<int>& dirichlet_colors() const { return colors; }
  virtual std::set<int>& dirichlet_colors() { return colors; }

  virtual const std::vector<int>& components_on_color(int c) const
  {
    assert(comp_on_color.find(c) != comp_on_color.end());
    return comp_on_color.find(c)->second;
  }
  virtual std::vector<int>& components_on_color(int c)
  {
    assert(comp_on_color.find(c) != comp_on_color.end());
    return comp_on_color.find(c)->second;
  }
};
} // namespace Gascoigne

#endif
