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

#include <array>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../Interface/gascoigne.h"

#include "nvector.h"

/***************************************************/

namespace Gascoigne {
class DataFormatHandler
{
private:
  typedef std::pair<std::string, std::string> NameType;
  typedef std::pair<std::string, double> StringDouble;

  template<typename T>
  using NameMap = std::map<std::string, T>;

  // typedef NameMap<void*> TypeVoid;
  typedef NameMap<std::string*> TypeString;
  typedef NameMap<int*> TypeInt;
  typedef NameMap<IndexType*> TypeIndex;
  typedef NameMap<bool*> TypeBool; // neu
  typedef NameMap<float*> TypeFloat;
  typedef NameMap<double*> TypeDouble;
  typedef NameMap<StringDouble*> TypeStringDouble;

  typedef NameMap<std::array<double, 2>*> TypeFix2Double;
  typedef NameMap<std::array<double, 3>*> TypeFix3Double;

  typedef NameMap<std::vector<double>*> TypeVectorDouble;
  typedef NameMap<IntVector*> TypeVectorInt;
  typedef NameMap<IndexVector*> TypeVectorIndex;
  typedef NameMap<std::vector<std::string>*> TypeVectorString;

  typedef NameMap<IntSet*> TypeSetInt;
  typedef NameMap<IndexSet*> TypeIndexSet;
  typedef NameMap<std::set<std::vector<std::string>>*> TypeSetVectorString;

  typedef NameMap<std::map<int, IntVector>*> TypeMapIntVectorInt;
  typedef NameMap<std::map<IndexType, IndexVector>*> TypeMapIndexVectorIndex;

  std::map<std::string, std::string> NT;

  // TypeVoid TV;

  TypeString TS;
  TypeInt TI;
  TypeIndex TId;
  TypeBool TB; // neu
  TypeFloat TF;
  TypeDouble TD;
  TypeFix2Double TF2D;
  TypeFix3Double TF3D;

  TypeVectorDouble TND;
  TypeVectorInt TNI;
  TypeVectorIndex TNId;
  TypeVectorString TVS;

  TypeSetInt TSI;
  TypeIndexSet TSId;
  TypeSetVectorString TSVS;

  TypeMapIntVectorInt TMINI;
  TypeMapIndexVectorIndex TMIdNId;

  TypeStringDouble TSD;

public:
  void clear()
  {
    NT.clear();
    TS.clear();
    TI.clear();
    TId.clear();
    TB.clear();
    TF.clear();
    TD.clear();
    TF2D.clear();
    TF3D.clear();
    TND.clear();
    TNI.clear();
    TNId.clear();
    TVS.clear();
    TSI.clear();
    TSId.clear();
    TSVS.clear();
    TMINI.clear();
    TMIdNId.clear();
    TSD.clear();
  }

  // template<typename T>
  // void insert(const std::string& name, T val);

  // template<typename T>
  // void insert(const std::string& name, T* val, T def);

  // template<typename T>
  // void setvalue(const std::string& name, T val);

  // without default values
  void insert(const std::string&, std::string*);
  void insert(const std::string&, int*);
  void insert(const std::string&, IndexType*);
  void insert(const std::string&, bool*); // neu
  void insert(const std::string&, float*);
  void insert(const std::string&, double*);

  void insert(const std::string&, std::array<double, 2>*);
  void insert(const std::string&, std::array<double, 3>*);

  void insert(const std::string&, std::vector<double>*);
  void insert(const std::string&, IntVector*);
  void insert(const std::string&, IndexVector*);
  void insert(const std::string&, std::vector<std::string>*);

  void insert(const std::string&, IntSet*);
  void insert(const std::string&, IndexSet*);
  void insert(const std::string&, std::set<std::vector<std::string>>*);

  void insert(const std::string&, std::map<int, IntVector>*);
  void insert(const std::string&, std::map<IndexType, IndexVector>*);
  void insert(const std::string&, std::map<int, std::string>*);

  void insert(int, StringDouble*);

  // with default values
  void insert(const std::string&, std::string*, const std::string&);
  void insert(const std::string&, int*, int);
  void insert(const std::string&, IndexType*, IndexType);
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
  void insert(const std::string&, IndexVector*, IndexVector&);

  void get(std::string&, const std::string&);

  void setvalue(const std::string&, const std::string&);
  void setvalue(const std::string&, int);
  void setvalue(const std::string&, IndexType);
  void setvalue(const std::string&, bool); // neu
  void setvalue(const std::string&, float);
  void setvalue(const std::string&, double);

  void setvalue(const std::string&, std::array<double, 2>&);
  void setvalue(const std::string&, std::array<double, 3>&);

  void setvalue(const std::string&, std::vector<double>&);
  void setvalue(const std::string&, IntVector&);
  void setvalue(const std::string&, IndexVector&);
  void setvalue(const std::string&, std::vector<std::string>&);

  void setvalue(const std::string&, IntSet&);
  void setvalue(const std::string&, IndexSet&);

  void setvalue(const std::string&, std::pair<int, IntVector>&);
  void setvalue(const std::string&, std::pair<IndexType, IndexVector>&);

  void setvalue(const std::string&, StringDouble&);

  void insertvalue(const std::string&, std::vector<std::string>&);

  void print(std::ostream&) const;
};
} // namespace Gascoigne

#endif
