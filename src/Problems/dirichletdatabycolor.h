/**
 *
 * Copyright (C) 2004, 2007, 2008 by the Gascoigne 3D authors
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

#ifndef __DirichletDataByColor_h
#define __DirichletDataByColor_h

#include "../Interface/dirichletdata.h"

/*-----------------------------------------*/

namespace Gascoigne {
class DirichletDataByColor : public DirichletData
{
protected:
  std::set<IndexType> __cols;
  nvector<IndexType> __comps;
  nvector<double> __scales;

public:
  DirichletDataByColor(const ParamFile& pf,
                       nvector<IndexType> comps,
                       const std::set<IndexType>& cl,
                       nvector<double> s);
  DirichletDataByColor(nvector<IndexType> comps,
                       const std::set<IndexType>& cl,
                       nvector<double> s);
  DirichletDataByColor(IndexType comps, std::set<IndexType>& cl, double s);
  DirichletDataByColor(const std::vector<std::string>& args);
  ~DirichletDataByColor() {}

  std::string GetName() const { return "DirichletDataByColor"; }

  std::set<IndexType> preferred_colors() const;
  void operator()(DoubleVector& b, const Vertex2d& V, IndexType color) const;
  void operator()(DoubleVector& b, const Vertex3d& V, IndexType color) const;
};
} // namespace Gascoigne

#endif
