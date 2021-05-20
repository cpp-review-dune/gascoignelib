
#ifndef __local_h
#define __local_h

#include "componentinformationbase.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "problemdescriptorbase.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;

template<int DIM>
class Fluid_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "Fluid CI"; }

  const int GetNScalars() const { return DIM + 2; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "p";
    if (DIM == 2) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
    }
    if (DIM == 3) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
      if (i == 3)
        s_name = "vz";
    }
  }

  const int GetNVectors() const { return 1; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "V";
    else
      abort();
  }
  void GetVectorIndices(int i, array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 1;
      fa_vectorindices[1] = 2;
      // fa_vectorindices[2]=-1;
      if (DIM == 3)
        fa_vectorindices[2] = 3;
    }
  }
};

template<int DIM>
class Solid_Euler_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "Solid_Euler CI"; }

  const int GetNScalars() const { return DIM + 1; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "p";
    if (DIM == 2) {
      if (i == 1)
        s_name = "ux";
      if (i == 2)
        s_name = "uy";
    }
    if (DIM == 3) {
      if (i == 1)
        s_name = "ux";
      if (i == 2)
        s_name = "uy";
      if (i == 3)
        s_name = "uz";
    }
  }

  const int GetNVectors() const { return 1; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "U";
    else
      abort();
  }
  void GetVectorIndices(int i, array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 1;
      fa_vectorindices[1] = 2;
      // fa_vectorindices[2]=-1;
      if (DIM == 3)
        fa_vectorindices[2] = 3;
    }
  }
};
template<int DIM>
class FSI_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "FSI CI"; }

  const int GetNScalars() const { return 2 * DIM + 2; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "p";
    if (DIM == 2) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
      if (i == 3)
        s_name = "ux";
      if (i == 4)
        s_name = "uy";
    }
    if (DIM == 3) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
      if (i == 3)
        s_name = "vz";
      if (i == 4)
        s_name = "ux";
      if (i == 5)
        s_name = "uy";
      if (i == 6)
        s_name = "uz";
    }
  }

  const int GetNVectors() const { return 2; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "V";
    else if (i == 1)
      s_name = "U";
    else
      abort();
  }
  void GetVectorIndices(int i, array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 1;
      fa_vectorindices[1] = 2;
      // fa_vectorindices[2]=-1;
      if (DIM == 3)
        fa_vectorindices[2] = 3;
    } else {
      fa_vectorindices[0] = DIM + 1;
      fa_vectorindices[1] = DIM + 2;
      // fa_vectorindices[2]=-1;
      if (DIM == 3)
        fa_vectorindices[2] = DIM + 3;
    }
  }
};

#endif
