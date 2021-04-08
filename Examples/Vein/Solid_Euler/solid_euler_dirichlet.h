
#ifndef __solid_euler_dirichlet_h
#define __solid_euler_dirichlet_h

#include "dirichletdata.h"
#include "filescanner.h"
#include "paramfile.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class DD_Solid_Euler2d : public DirichletData {
protected:
public:
  DD_Solid_Euler2d(const ParamFile *pf) {}

  std::string GetName() const { return "DD_Solid_Euler"; }

  void operator()(DoubleVector &b, const Vertex2d &v, int color) const {
    b.zero();
  }
};

class DD_Solid_Euler3d : public DirichletData {
protected:
public:
  DD_Solid_Euler3d(const ParamFile *pf) {}

  std::string GetName() const { return "DD_Solid_Euler3d"; }

  void operator()(DoubleVector &b, const Vertex3d &v, int color) const {
    b.zero();
  }
};

#endif
