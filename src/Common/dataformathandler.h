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

#ifndef __dataformathandler_h
#define __dataformathandler_h

#include "nvector.h"
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "gascoigne.h"
#include <array>
#include <string>

/***************************************************/

namespace Gascoigne {
class DataFormatHandler
{
private:
  typedef std::pair<std::string, std::string> NameType;
  typedef std::pair<std::string, double> StringDouble;

  template<typename T>
  using NameMap = std::map<std::string, T>;

  typedef NameMap<std::string*> TypeString;
  typedef NameMap<int*> TypeInt;
  typedef NameMap<bool*> TypeBool; // neu
  typedef NameMap<float*> TypeFloat;
  typedef NameMap<double*> TypeDouble;
  typedef NameMap<StringDouble*> TypeStringDouble;

  typedef NameMap<std::array<double, 2>*> TypeFix2Double;
  typedef NameMap<std::array<double, 3>*> TypeFix3Double;

  typedef NameMap<std::vector<double>*> TypeVectorDouble;
  typedef NameMap<IntVector*> TypeVectorInt;
  typedef NameMap<std::vector<std::string>*> TypeVectorString;

  typedef NameMap<std::set<int>*> TypeSetInt;
  typedef NameMap<std::set<std::vector<std::string>>*> TypeSetVectorString;

  typedef NameMap<std::map<int, IntVector>*> TypeMapIntVectorInt;

  std::set<NameType> NT;

  TypeString TS;
  TypeInt TI;
  TypeBool TB; // neu
  TypeFloat TF;
  TypeDouble TD;
  TypeFix2Double TF2D;
  TypeFix3Double TF3D;

  TypeVectorDouble TND;
  TypeVectorInt TNI;
  TypeVectorString TVS;

  TypeSetInt TSI;
  TypeSetVectorString TSVS;

  TypeMapIntVectorInt TMINI;

  TypeStringDouble TSD;

  std::string search(std::string& fo, const std::string& name);

public:
  void clear()
  {
    NT.clear();
    TS.clear();
    TI.clear();
    TB.clear();
    TF.clear();
    TD.clear();
    TF2D.clear();
    TF3D.clear();
    TND.clear();
    TNI.clear();
    TVS.clear();
    TSI.clear();
    TSVS.clear();
    TMINI.clear();
    TSD.clear();
  }

  // without default values
  void insert(const std::string&, std::string*);
  void insert(const std::string&, int*);
  void insert(const std::string&, bool*); // neu
  void insert(const std::string&, float*);
  void insert(const std::string&, double*);

  void insert(const std::string&, std::array<double, 2>*);
  void insert(const std::string&, std::array<double, 3>*);

  void insert(const std::string&, std::vector<double>*);
  void insert(const std::string&, IntVector*);
  void insert(const std::string&, std::vector<std::string>*);

  void insert(const std::string&, IntSet*);
  void insert(const std::string&, std::set<std::vector<std::string>>*);

  void insert(const std::string&, std::map<int, IntVector>*);
  void insert(const std::string&, std::map<int, std::string>*);

  void insert(int, StringDouble*);

  // with default values
  void insert(const std::string&, std::string*, const std::string&);
  void insert(const std::string&, int*, int);
  void insert(const std::string&, bool*, bool); // neu
  void insert(const std::string&, float*, float);
  void insert(const std::string&, double*, double);
  void insert(const std::string&,
              std::array<double, 2>*,
              std::array<double, 2>&);
  void insert(const std::string&,
              std::array<double, 3>*,
              std::array<double, 3>&);
  void insert(const std::string&, std::vector<double>*, std::vector<double>&);
  void insert(const std::string&, IntVector*, IntVector&);

  void get(std::string&, const std::string&);

  void setvalue(const std::string&, const std::string&);
  void setvalue(const std::string&, int);
  void setvalue(const std::string&, bool); // neu
  void setvalue(const std::string&, float);
  void setvalue(const std::string&, double);

  void setvalue(const std::string&, std::array<double, 2>&);
  void setvalue(const std::string&, std::array<double, 3>&);

  void setvalue(const std::string&, std::vector<double>&);
  void setvalue(const std::string&, IntVector&);
  void setvalue(const std::string&, std::vector<std::string>&);

  void setvalue(const std::string&, IntSet&);

  void setvalue(const std::string&, std::pair<int, IntVector>&);

  void setvalue(const std::string&, StringDouble&);

  void insertvalue(const std::string&, std::vector<std::string>&);

  void print(std::ostream&) const;
};
} // namespace Gascoigne

#endif
