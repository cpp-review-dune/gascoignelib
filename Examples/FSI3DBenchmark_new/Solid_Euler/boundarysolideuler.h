/*----------------------------    boundarysolideuler.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __BoundarySolidEuler_H
#define __BoundarySolidEuler_H
/*----------------------------    boundarysolideuler.h
 * ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> class BoundarySolidEuler : public BoundaryEquation {

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  mutable VECTOR __n, PHI;
  mutable MATRIX F, NU, NPHI;
  mutable double J, __pressure;
  // mutable FemFunction *OLD, *DEF, *DEFOLD;

  mutable MATRIX NV;
  mutable VECTOR g, V;
  double __nu_f, __rho_f;
  double p_2, p_4;
  mutable FemFunction *VEL;

  void SetFemData(FemData &q) const {
    if (q.find("VEL") == q.end()) {
      std::cout << "VEL in solideulerbound not found" << std::endl;
      abort();
    }
    VEL = &q["VEL"];
  }

public:
  ~BoundarySolidEuler() {}
  BoundarySolidEuler() { abort(); }
  BoundarySolidEuler(const ParamFile *pf);

  std::string GetName() const { return "BoundaryPrestress"; }

  int GetNcomp() const { return DIM; }
  void Form(VectorIterator b, const FemFunction &U, const TestFunction &N,
            int col) const;
  void Matrix(EntryMatrix &E, const FemFunction &U, const TestFunction &M,
              const TestFunction &N, int col) const;
  void pointboundary(double h, const FemFunction &U, const Vertex<DIM> &v,
                     const Vertex<DIM> &n) const;
};

} // namespace Gascoigne

/*----------------------------   boundarysolideuler.h
 * ---------------------------*/
/* end of #ifndef __ boundarysolideuler_H */
#endif
/*----------------------------    boundarysolideuler.h
 * ---------------------------*/
