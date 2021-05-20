/*----------------------------    boundaryfsi.h ---------------------------*/
/*      $Id:$                 */
#ifndef __BoundaryFSI_H
#define __BoundaryFSI_H
/*----------------------------    boundaryfsi.h ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template<int DIM>
class BoundaryFSI : public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  mutable VECTOR g, g_OLD;
  mutable VECTOR __n, PHI;
  mutable MATRIX F, NV, NU, F_OLD, NV_OLD, NU_OLD, NPHI;
  mutable double J, J_OLD;
  double __nu_f, __rho_f;
  double p_2, p_4;

  mutable FemFunction *U_Vec, *UOLD_Vec;

  void SetFemData(FemData& q) const
  {

    assert(q.find("U_Vec") != q.end());
    U_Vec = &q["U_Vec"];

    assert(q.find("UOLD_Vec") != q.end());
    UOLD_Vec = &q["UOLD_Vec"];
  }

public:
  ~BoundaryFSI() {}
  BoundaryFSI() { abort(); }
  BoundaryFSI(const ParamFile* pf);

  std::string GetName() const { return "BoundaryFSI"; }

  int GetNcomp() const { return DIM + 1; }

  void Form(VectorIterator b,
            const FemFunction& U,
            const TestFunction& N,
            int col) const;
  void Matrix(EntryMatrix& E,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N,
              int col) const;

  void pointboundary(double h,
                     const FemFunction& U,
                     const Vertex<DIM>& v,
                     const Vertex<DIM>& n) const;
};

} // namespace Gascoigne

/*----------------------------   boundaryfsi.h     ---------------------------*/
/* end of #ifndef __ boundaryfsi_H */
#endif
/*----------------------------    boundaryfsi.h ---------------------------*/
