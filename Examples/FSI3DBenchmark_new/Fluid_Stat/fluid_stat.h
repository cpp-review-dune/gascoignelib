/*----------------------------   fluid_stat.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fluid_stat_H
#define __fluid_stat_H
/*----------------------------   fluid_stat.h     ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"
#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM>
class Fluid_Stat : public LpsEquation // , public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  // problem parameter
  double rho_f, nu_f;
  double lps0;

  mutable double __h, domain;
  mutable Vertex<DIM> __v;

  // stuff from point_M
  mutable double divergence;
  mutable MATRIX CONV_dV1, PRESSURE_P;
  mutable std::array<double, DIM> DOMAIN_V;
  mutable std::array<double, DIM> DIVERGENCE_V;
  mutable std::array<double, DIM> CONV_dV2;

  mutable MATRIX TENSOR_dV[DIM];
  // stuff from point
  mutable double lps;
  mutable MATRIX NV;
  mutable MATRIX SIGMAf;

  mutable VECTOR V;

public:
  ~Fluid_Stat() {}
  Fluid_Stat() { abort(); }
  Fluid_Stat(const ParamFile *pf);

  std::string GetName() const { return "Fluid_Stat"; }

  int GetNcomp() const { return DIM + 1; }

  void point(double h, const FemFunction &U, const Vertex<DIM> &v) const;
  void point_M(int j, const FemFunction &U, const TestFunction &M) const;
  void point_cell(int material) const;

  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const;

  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const;
  void MatrixBlock(EntryMatrix &A, const FemFunction &U,
                   const FemFunction &NNN) const;

  ////////////////////////////////////////////////// LPS

  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const;

  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &N) const;

  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const;
};

} // namespace Gascoigne

/*----------------------------   fluid_stat.h     ---------------------------*/
/* end of #ifndef __fluid_stat_H */
#endif
/*----------------------------   fluid_stat.h     ---------------------------*/
