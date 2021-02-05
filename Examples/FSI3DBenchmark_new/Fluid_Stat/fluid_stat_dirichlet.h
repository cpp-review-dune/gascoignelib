
#ifndef __fluid_stat_dirichlet_h
#define __fluid_stat_dirichlet_h

#include "dirichletdata.h"
#include "filescanner.h"
#include "paramfile.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

class DD_Fluid_Stat2d : public DirichletData {
protected:
  double vmean;

public:
  DD_Fluid_Stat2d(const ParamFile *pf) {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector &b, const Vertex2d &v, int color) const {
    b.zero();
  }
};

class DD_Fluid_Stat3d : public DirichletData {
protected:
  double vmean;

public:
  DD_Fluid_Stat3d(const ParamFile *pf) {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
    cout << "%%%%%%%%%%Fluid_Stat%%%%%%%%%%" << endl;
    cout << "Boundary 7 -- inflow 100ml/min" << endl;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  }

  std::string GetName() const { return "DD_Fluid_Stat3d"; }

  void operator()(DoubleVector &b, const Vertex3d &v, int color) const {
    b.zero();
    // inflow radius
    double radius = 0.3;
    // inflow in ml/s
    double flow = 100. / 60.;
    double fact_time = 2 * flow / (radius * radius * M_PI);

    double x = v.x();
    double y = v.y();
    double z = v.z();

    if (color == 7) {
      Vertex3d center;
      center[0] = 3.78223304703;
      center[1] = -3.882233047;
      center[2] = 1.0928932;
      Vertex3d normal;
      normal[0] = 0.5000;
      normal[1] = -0.4998;
      normal[2] = 0.7072;

      b[1] += -normal[0] * fact_time *
              ((x - center[0]) * (x - center[0]) +
               (y - center[1]) * (y - center[1]) +
               (z - center[2]) * (z - center[2]) - radius * radius) /
              (-radius * radius);
      b[2] += -normal[1] * fact_time *
              ((x - center[0]) * (x - center[0]) +
               (y - center[1]) * (y - center[1]) +
               (z - center[2]) * (z - center[2]) - radius * radius) /
              (-radius * radius);
      b[3] += -normal[2] * fact_time *
              ((x - center[0]) * (x - center[0]) +
               (y - center[1]) * (y - center[1]) +
               (z - center[2]) * (z - center[2]) - radius * radius) /
              (-radius * radius);
    }
  }
};

#endif
