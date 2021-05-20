/*----------------------------   prestress.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solideuler_H
#define __solideuler_H
/*----------------------------   prestress.h     ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template<int DIM>
class SolidEuler : public LpsEquation // , public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  // problem parameter
  double rho_s, lambda_s, mu_s, kapa_s, rho_f, nu_f;
  std::string mat_law;
  // solid
  mutable MATRIX sigmas_dU[DIM];
  mutable double domain;
  // stuff from point
  mutable double __j;
  mutable MATRIX __e, SIGMAs, sigmas, __f, SIGMAs_dU, __c, NV, SIGMAf;
  mutable VECTOR V;
  double lps0;
  mutable double lps;
  // mutable FemFunction *OLD, *DEF, *DEFOLD;
  mutable FemFunction* VEL;

  void SetFemData(FemData& q) const
  {
    if (q.find("VEL") == q.end()) {
      std::cout << "VEL in solideuler not found" << std::endl;
      abort();
    }
    VEL = &q["VEL"];

    /*	assert(q.find("OLD")!=q.end());
            OLD = &q["OLD"];

            assert(q.find("DEF")!=q.end());
            DEF = &q["DEF"];

            assert(q.find("DEFOLD")!=q.end());
            DEFOLD = &q["DEFOLD"];
            */
  }

public:
  ~SolidEuler() {}
  SolidEuler() { abort(); }
  SolidEuler(const ParamFile* pf);

  std::string GetName() const { return "Prestress"; }

  int GetNcomp() const { return DIM + 1; }

  void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;
  void point_M(int j, const FemFunction& U, const TestFunction& M) const;

  void Form(VectorIterator b,
            const FemFunction& U,
            const TestFunction& N) const;

  void Matrix(EntryMatrix& A,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N) const;
  void MatrixBlock(EntryMatrix& A,
                   const FemFunction& U,
                   const FemFunction& NNN) const;
  void point_cell(int material) const;

  ////////////////////////////////////////////////// LPS

  void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const {}

  void StabForm(VectorIterator b,
                const FemFunction& U,
                const FemFunction& UP,
                const TestFunction& N) const
  {}

  void StabMatrix(EntryMatrix& A,
                  const FemFunction& U,
                  const TestFunction& Np,
                  const TestFunction& Mp) const
  {}
};

} // namespace Gascoigne

/*----------------------------   solideuler.h     ---------------------------*/
/* end of #ifndef __SolidEuler_H */
#endif
/*----------------------------   solideuler.h     ---------------------------*/
