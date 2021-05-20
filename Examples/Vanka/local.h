
#ifndef __local_h
#define __local_h

#include "componentinformationbase.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "fsi.h"
#include "problemdescriptorbase.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;

template<int DIM>
class FSI_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "FSI CI"; }

  const int GetNScalars() const { return DIM; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "ux";
    if (i == 1)
      s_name = "uy";
    if (DIM == 3)
      if (i == 2)
        s_name = "uz";
  }

  const int GetNVectors() const { return 1; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "U";
    else
      abort();
  }
  void GetVectorIndices(int i, std::array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 0;
      fa_vectorindices[1] = 1;
      fa_vectorindices[2] = -1;
      if (DIM == 3)
        fa_vectorindices[2] = 2;
    }
  }
};

template<int DIM>
class MyDD : public DirichletData
{
protected:
  double vmean;

public:
  MyDD(const ParamFile* pf) {}

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();
  }
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
  }
};

template<int DIM>
class ProblemDescriptor : public ProblemDescriptorBase
{
public:
  std::string GetName() const { return "fsi"; }
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new Gascoigne::EQ<DIM>(GetParamFile());
    GetDirichletDataPointer() = new MyDD<DIM>(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    GetComponentInformationPointer() = new FSI_CI<DIM>;
  }
};

#endif
