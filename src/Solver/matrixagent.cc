/**
 *
 * Copyright (C) 2004, 2005, 2006, 2020 by the Gascoigne 3D authors
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

#include "matrixagent.h"

namespace Gascoigne {

/*-------------------------------------------------*/

MatrixAgent::MatrixAgent() {}

/*-------------------------------------------------*/

MatrixAgent::~MatrixAgent()
{
  for (auto p = begin(); p != end(); p++) {
    if (p->second) {
      delete p->second;
      p->second = NULL;
    }
  }
  std::map<Matrix, MatrixInterface*>::clear();
}

/*-------------------------------------------------*/

void
MatrixAgent::Register(const Matrix& mg)
{
  auto p = find(mg);
  if (p == end())
    insert(std::make_pair(mg, static_cast<MatrixInterface*>(NULL)));
}

/*-------------------------------------------------*/

void
MatrixAgent::Delete(Matrix& mg)
{
  auto p = find(mg);
  if (p != end()) {
    delete p->second;
    erase(p);
  }
}

/*-------------------------------------------------*/

MatrixInterface&
MatrixAgent::operator()(const Matrix& g)
{
  auto p = find(g);
  if (p == end()) {
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": MatrixAgent::operator(): ERROR" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": Matrix '" << g << "' not found in list of: " << std::endl;
    std::cerr << " " << *this << std::endl;
    abort();
  }
  MatrixInterface* vp = p->second;
  if (vp == NULL) {
    std::cerr << "MatrixAgent  MatrixInterface* NULL\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    abort();
  }
  return *vp;
}

const MatrixInterface&
MatrixAgent::operator()(const Matrix& g) const
{
  const auto p = find(g);
  if (p == end()) {
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": MatrixAgent::operator(): ERROR" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": Matrix '" << g << "' not found in list of: " << std::endl;
    std::cerr << " " << *this << std::endl;
    abort();
  }
  MatrixInterface* vp = p->second;
  if (vp == NULL) {
    std::cerr << "MatrixAgent  MatrixInterface* NULL\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    abort();
  }
  return *vp;
}

} // namespace Gascoigne
