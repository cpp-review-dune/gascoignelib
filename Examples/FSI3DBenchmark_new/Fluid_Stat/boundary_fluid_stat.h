/*----------------------------    boundary_fluid_stat.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __Bboundary_fluid_stat_H
#define __boundary_fluid_stat_H
/*----------------------------    boundary_fluid_stat.h
 * ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"
#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> class Boundary_Fluid_Stat : public BoundaryEquation {

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  mutable VECTOR g;
  mutable VECTOR __n, PHI;
  mutable MATRIX NV, NPHI;

  double __nu_f, __rho_f;
  double p_2, p_4;

public:
  ~Boundary_Fluid_Stat() {}
  Boundary_Fluid_Stat() { abort(); }
  Boundary_Fluid_Stat(const ParamFile *pf);

  std::string GetName() const { return "Boundary_Fluid_Stat"; }

  int GetNcomp() const { return DIM + 1; }

  void Form(VectorIterator b, const FemFunction &U, const TestFunction &N,
            int col) const;
  void Matrix(EntryMatrix &E, const FemFunction &U, const TestFunction &M,
              const TestFunction &N, int col) const;

  void pointboundary(double h, const FemFunction &U, const Vertex<DIM> &v,
                     const Vertex<DIM> &n) const;
};

} // namespace Gascoigne

/*----------------------------   boundary_fluid_stat.h
 * ---------------------------*/
/* end of #ifndef __boundary_fluid_stat_H */
#endif
/*----------------------------    boundary_fluid_stat.h
 * ---------------------------*/
