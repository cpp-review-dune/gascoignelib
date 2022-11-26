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

#ifndef __DiracRightHandSide_h
#define __DiracRightHandSide_h

#include <vector>

#include "../Common/compvector.h"
#include "../Common/vertex.h"

#include "application.h"

/**********************************************************/

namespace Gascoigne {
class DiracRightHandSide : public virtual Application
{
private:
protected:
  mutable std::vector<Vertex2d> _v2d;
  mutable std::vector<Vertex3d> _v3d;
  mutable std::vector<ShortIndexType> _comps;
  mutable std::vector<double> _weights;

public:
  DiracRightHandSide() {}
  virtual ~DiracRightHandSide() {}

  /**
     clones a DRHS. Usually it is simply to return a new instance of the same
     object passing the required variables to the constructor. It takes the role
     of a copy constructor and the cloning of classes is required for
     multithreading.
  */
  virtual DiracRightHandSide* createNew() const
  {
    std::cerr << "\"DRHS::createNew\" not written!" << std::endl;
    abort();
  }

  virtual void BasicInit(const std::vector<Vertex2d>& v2d,
                         const std::vector<ShortIndexType>& comps,
                         const std::vector<double>& weights)
  {
    _v2d = v2d;
    _comps = comps;
    _weights = weights;
  }
  virtual void BasicInit(const std::vector<Vertex3d>& v3d,
                         const std::vector<ShortIndexType>& comps,
                         const std::vector<double>& weights)
  {
    _v3d = v3d;
    _comps = comps;
    _weights = weights;
  }

  virtual const std::vector<Vertex2d>& GetPoints2d() const { return _v2d; }
  virtual const std::vector<Vertex3d>& GetPoints3d() const { return _v3d; }

  virtual const std::vector<ShortIndexType>& GetComps() const { return _comps; }

  virtual const std::vector<double>& GetWeights() const { return _weights; }

  virtual double operator()(int i, const Vertex2d& v) const
  {
    std::cerr << "\"DiracRightHandSide::operator()\" not written!" << std::endl;
    abort();
  }
  virtual double operator()(int i, const Vertex3d& v) const
  {
    std::cerr << "\"DiracRightHandSide::operator()\" not written!" << std::endl;
    abort();
  }

  virtual void operator()(int i,
                          VectorIterator b,
                          const TestFunction& N,
                          const Vertex2d& v) const
  {
    b[_comps[i]] += N.m() * (*this)(i, v);
  }
  virtual void operator()(int i,
                          VectorIterator b,
                          const TestFunction& N,
                          const Vertex3d& v) const
  {
    b[_comps[i]] += N.m() * (*this)(i, v);
  }
};

typedef DiracRightHandSide DiracInitialCondition;

/**********************************************************/

} // namespace Gascoigne

#endif
