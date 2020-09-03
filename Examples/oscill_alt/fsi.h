/*----------------------------   fsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsi_H
#define __fsi_H
/*----------------------------   fsi.h     ---------------------------*/
#define BDF
#include "equation.h"
#include "boundaryequation.h"
#include "paramfile.h"
#include "lpsequation.h"
#include <eigen3/Eigen/Dense>
#include "chi.h"

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM>
class FSI : public LpsEquation  //, public BoundaryEquation
{
protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1>   VECTOR;

  // problem parameter TODO add __Theta and setTheta method
  double         rho_f, nu_f, rho_s, lambda_s, mu_s, nu_e;
  mutable double theta;
  double         extend0, lps0;
  double         pp_e;

  mutable double      __h;
  mutable Vertex<DIM> __v;

  // stuff from point_M
  mutable double                  divergence;
  mutable MATRIX                  CONV_dV1, DOMAIN_U1, PRESSURE_P;
  mutable VECTOR                  DOMAIN_U2;
  mutable std::array<double, DIM> DOMAIN_V;
  mutable std::array<double, DIM> Jj, DIVERGENCE_U, DIVERGENCE_V;
  mutable std::array<double, DIM> CONV_dV2;
  mutable MATRIX                  Fij[DIM];
  mutable MATRIX                  TENSOR_dU[DIM];
  mutable MATRIX                  TENSOR_dV[DIM];
  mutable MATRIX                  PRESSURE_U[DIM];
  mutable VECTOR                  CONV_dU[DIM];

  mutable MATRIX DAMP, DAMP_old;
  double         damp0;

  // solid
  mutable VECTOR SIGMA_dF[DIM];
  mutable MATRIX SIGMA_dU[DIM];

  // boundary
  // mutable VECTOR normal;
  // mutable VECTOR BOUNDARY, BOUNDARY_old;

  // stuff from point
  mutable double extend, lps, J, J_old;
  mutable MATRIX NV, NU, NV_old, NU_old, F, F_old, SIGMAs, SIGMAs_old, SIGMAs_jU, E;
  mutable MATRIX SIGMAf, SIGMAf_old;

  mutable VECTOR V, V_old, dtU;

  Chi                  chi;
  mutable int          domain;
  mutable bool         bdf1_= true;
  mutable FemFunction* OLD;
  mutable FemFunction* OLDOLD;
  void                 SetFemData(FemData& q) const override final {
    auto qold= q.find("OLD");
    if (qold != q.end()) {
      OLD= &qold->second;
    }
    qold= q.find("OLD0");
    if (qold != q.end()) {
      OLDOLD= &qold->second;
    }
  }

public:
  ~FSI() {
  }
  FSI() {
    abort();
  }
  FSI(const ParamFile* pf);

  std::string GetName() const override final {
    return "FSI";
  }

  int GetNcomp() const override final {
    return 2 * DIM + 1;
  }
  double& getTheta() const {
    return this->theta;
  }
  bool& bdf1() const {
    return this->bdf1_;
  }

  void point(double h, const FemFunction& U, const Vertex<DIM>& v) const override final;
  void point_M(int j, const FemFunction& U, const TestFunction& M) const;

  void Form(VectorIterator      b,
            const FemFunction&  U,
            const TestFunction& N) const override final;

  void Matrix(EntryMatrix&        A,
              const FemFunction&  U,
              const TestFunction& M,
              const TestFunction& N) const override final;
  void MatrixBlock(EntryMatrix&       A,
                   const FemFunction& U,
                   const FemFunction& NNN) const override final;

  /* ////////////////////////////////////////////////// Boundary */

  // void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col)
  // const; void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const
  // TestFunction& N,
  //             int col) const;
  //
  // void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v,
  //                    const Vertex<DIM>& n) const;

  ////////////////////////////////////////////////// LPS

  void
    lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const override final;

  void StabForm(VectorIterator      b,
                const FemFunction&  U,
                const FemFunction&  UP,
                const TestFunction& N) const override final;

  void StabMatrix(EntryMatrix&        A,
                  const FemFunction&  U,
                  const TestFunction& Np,
                  const TestFunction& Mp) const override final;
};

}  // namespace Gascoigne

/*----------------------------   fsi.h     ---------------------------*/
/* end of #ifndef __fsi_H */
#endif
/*----------------------------   fsi.h     ---------------------------*/
