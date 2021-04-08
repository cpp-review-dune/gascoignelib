#include "solid.h"
#include "filescanner.h"

#include "multiplex.xx"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> Solid<DIM>::Solid(const ParamFile *pf) {
  DataFormatHandler DFH;
  DFH.insert("rho_s", &rho_s, 0.0);
  DFH.insert("lambda_s", &lambda_s, 0.0);
  DFH.insert("mu_s", &mu_s, 0.0);
  DFH.insert("mat_law", &mat_law, "STVK");

  FileScanner FS(DFH, pf, "Equation");
  assert(rho_s > 0);
  assert(lambda_s > 0);
  assert(mu_s > 0);

  kapa_s = lambda_s + 2.0 / 3.0 * mu_s;
  cout << "%%%%%%%%%% Problem %%%%%%%%%%" << endl
       << "  rho_s / mu_s / lambda_s: " << rho_s << " / " << mu_s << " / "
       << lambda_s << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
}

//////////////////////////////////////////////////

template <int DIM>
void Solid<DIM>::point(double h, const FemFunction &U,
                       const Vertex<DIM> &v) const {
  //__h = h;
  //__v = v;
  multiplex_init_F<DIM>(F, U);

  J = F.determinant();
  if (mat_law == "STVK") {
    E = 0.5 * (F.transpose() * F - MATRIX::Identity());
    SIGMAs = (2.0 * mu_s * E + lambda_s * E.trace() * MATRIX::Identity());
  } else if (mat_law == "artery") {
    C = F.transpose() * F;
    SIGMAs = mu_s * pow(J, -2.0 / 3.0) *
                 (MATRIX::Identity() - 1.0 / 3.0 * C.trace() * C.inverse()) +
             0.5 * kapa_s * (pow(J, 2.0) - 1.0) * C.inverse();
  } else
    abort();
}

template <int DIM>
void Solid<DIM>::Form(VectorIterator b, const FemFunction &U,
                      const TestFunction &N) const {
  // phi =nabla N;
  VECTOR phi;
  multiplex_init_test<DIM>(phi, N);

  // FULL tensor F Sigma
  for (int i = 0; i < DIM; ++i)
    b[i] += (F * SIGMAs * phi)(i, 0);
}

template <int DIM>
void Solid<DIM>::Matrix(EntryMatrix &A, const FemFunction &U,
                        const TestFunction &M, const TestFunction &N) const {

  VECTOR phi;
  multiplex_init_test<DIM>(phi, N);
  VECTOR psi;
  multiplex_init_test<DIM>(psi, M);

  // F Sigma
  // wrt F,
  for (int j = 0; j < DIM; ++j)
    A(j, j) += (SIGMA_dF[j].transpose() * phi)(0, 0);

  for (int j = 0; j < DIM; ++j) // F Sigma_j \nabla phi
  {
    VECTOR X = SIGMA_dU[j] * phi;
    for (int i = 0; i < DIM; ++i) {
      A(i, j) += X(i, 0);
    }
  }
}

template <int DIM>
void Solid<DIM>::point_M(int j, const FemFunction &U,
                         const TestFunction &M) const {
  VECTOR psi;
  multiplex_init_test<DIM>(psi, M);
  if (mat_law == "STVK") {
    for (int j = 0; j < DIM; ++j)
      SIGMA_dF[j] = (psi.transpose() * SIGMAs).transpose();
    for (int j = 0; j < DIM; ++j) {
      MATRIX Ej = 0.5 * (psi * (F.block(j, 0, 1, DIM)) +
                         F.block(j, 0, 1, DIM).transpose() * psi.transpose());
      SIGMA_dU[j] =
          F * (2.0 * mu_s * Ej + lambda_s * Ej.trace() * MATRIX::Identity());
    }
  } else if (mat_law == "artery") {
    for (int j = 0; j < DIM; ++j)
      SIGMA_dF[j] = (psi.transpose() * SIGMAs).transpose();
    for (int j = 0; j < DIM; ++j) {
      MATRIX C_dU = (psi * (F.block(j, 0, 1, DIM)) +
                     F.block(j, 0, 1, DIM).transpose() * psi.transpose());
      MATRIX C_dU_inverse = -C.inverse() * C_dU * C.inverse();
      double J_dU =
          J * (F.inverse().block(0, j, DIM, 1) * psi.transpose()).trace();

      SIGMA_dU[j] = F * mu_s * (-2.0 / 3.0) * pow(J, -5.0 / 3.0) * J_dU *
                    (MATRIX::Identity() - 1.0 / 3.0 * C.trace() * C.inverse());
      SIGMA_dU[j] += F * mu_s * pow(J, -2.0 / 3.0) *
                     (-1.0 / 3.0 * C_dU.trace() * C.inverse());
      SIGMA_dU[j] += F * mu_s * pow(J, -2.0 / 3.0) *
                     (-1.0 / 3.0 * C.trace() * C_dU_inverse);
      SIGMA_dU[j] += F * 0.5 * kapa_s * 2.0 * J * J_dU * C.inverse();
      SIGMA_dU[j] += F * 0.5 * kapa_s * (pow(J, 2.0) - 1.0) * C_dU_inverse;
    }
  } else
    abort();
}

template <int DIM>
void Solid<DIM>::MatrixBlock(EntryMatrix &A, const FemFunction &U,
                             const FemFunction &N) const {
  ;
  for (int j = 0; j < N.size(); ++j) // trial
  {
#define M N[j]
    point_M(j, U, M);

    for (int i = 0; i < N.size(); ++i) // test
    {
      A.SetDofIndex(i, j);
      Matrix(A, U, M, N[i]);
    }
#undef M
  }
}

template class Solid<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
