/**
 *
 * Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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

#include "ghostagent.h"
#include "stlio.h"

namespace Gascoigne {

/*-------------------------------------------------*/

template<typename T>
GhostAgent<T>::GhostAgent()
{}

/*-------------------------------------------------*/

template<typename T>
GhostAgent<T>::~GhostAgent()
{
  for (auto p = std::map<std::string, T*>::begin();
       p != std::map<std::string, T*>::end();
       p++) {
    if (p->second) {
      delete p->second;
      p->second = NULL;
    }
  }
}

/*-------------------------------------------------*/

template<typename T>
void
GhostAgent<T>::Register(const std::string& mg)
{
  auto p = std::map<std::string, T*>::find(mg);
  if (p == std::map<std::string, T*>::end()) {
    std::map<std::string, T*>::emplace(mg, static_cast<T*>(NULL));
  }
}

/*-------------------------------------------------*/

template<typename T>
void
GhostAgent<T>::Delete(std::string& mg)
{
  auto p = std::map<std::string, T*>::find(mg);
  if (p != std::map<std::string, T*>::end()) {
    delete p->second;
    std::map<std::string, T*>::erase(p);
  }
}

/*-------------------------------------------------*/

template<typename T>
T&
GhostAgent<T>::operator()(const std::string& g)
{
  auto p = std::map<std::string, T*>::find(g);
  if (p == std::map<std::string, T*>::end()) {
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": GhostAgent::operator(): ERROR" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": Ghostvector '" << g
              << "' not found in list of: " << std::endl;
    std::cerr << " " << *this << std::endl;
    abort();
  }
  T* vp = p->second;
  if (vp == NULL) {
    std::cerr << "GhostAgent  T* NULL\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    abort();
  }
  return *vp;
}

template<typename T>
const T&
GhostAgent<T>::operator()(const std::string& g) const
{
  auto p = std::map<std::string, T*>::find(g);
  if (p == std::map<std::string, T*>::end()) {
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": GhostAgent::operator(): ERROR" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__;
    std::cerr << ": Ghostvector '" << g
              << "' not found in list of: " << std::endl;
    std::cerr << " " << *this << std::endl;
    abort();
  }
  T* vp = p->second;
  if (vp == NULL) {
    std::cerr << "GhostAgent  T* NULL\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    abort();
  }
  return *vp;
}

template class GhostAgent<GlobalVector>;
template class GhostAgent<MatrixInterface>;

} // namespace Gascoigne
