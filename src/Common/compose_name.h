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

#ifndef __compose_name_h
#define __compose_name_h

#include <string>

#include "../Interface/gascoigne.h"

namespace Gascoigne {
void
compose_name(std::string&, double, std::string f = "%5.2f");
void
compose_name(std::string&, IndexType);
void
compose_name(std::string&, IndexType, std::string t);
void
compose_name(std::string&, IndexType, IndexType);
void
compose_name_without_dot(std::string&, IndexType);
} // namespace Gascoigne

#endif
