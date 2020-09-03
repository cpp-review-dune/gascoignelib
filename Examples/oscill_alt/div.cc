#include "div.h"
#include "filescanner.h"
#include "multiplex.h"
extern double __DT;
extern double __THETA;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM>
DIV_proj<DIM>::DIV_proj(const ParamFile* pf) {
  DataFormatHandler DFH;
  DFH.insert("rho_s", &rho_s, 0.0);
  DFH.insert("pp_e", &pp_e, 0.0);
  DFH.insert("lambda_s", &lambda_s, 0.0);
  DFH.insert("mu_s", &mu_s, 0.0);
  DFH.insert("rho_f", &rho_f, 0.0);
  DFH.insert("nu_f", &nu_f, 0.0);
  DFH.insert("extend", &extend0, 0.0);
  DFH.insert("lps", &lps0, 0.0);
  DFH.insert("nu_e", &nu_e, 0.0);

  FileScanner FS(DFH, pf, "Equation");
  assert(rho_f > 0);
  assert(nu_f > 0);
  assert(extend0 > 0);
  assert(lps0 > 0);
  assert(rho_s > 0);
  assert(lambda_s > 0);
  assert(mu_s > 0);

  cout << "%%%%%%%%%% Problem %%%%%%%%%%" << endl
       << "  rho_f / nu_f: " << rho_f << " / " << nu_f << endl
       << "  rho_s / mu_s / lambda_s: " << rho_s << " / " << mu_s << " / " << lambda_s
       << endl
       << "  extend / lps: " << extend0 << " / " << lps0 << endl;

  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
}

//////////////////////////////////////////////////

template <int DIM>
void DIV_proj<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const {
  // domain = -1;
  domain= chi(v);
  if (domain < 0) {
    __h= h;
    __v= v;

    // set F, F_old to Identity in fluid
    Multiplex::init_F<DIM>(F, *OLD);
    Multiplex::init_F<DIM>(F_old, *OLD);

    auto F_inv    = F.inverse();
    auto F_inv_old= F_old.inverse();

    J    = F.determinant();
    J_old= F_old.determinant();

    Multiplex::init_NV<DIM>(NV, U);

    Multiplex::init_NV<DIM>(NV_old, *OLD);

    // TENSOR
    SIGMAf= rho_f * nu_f * (F_inv * NV + NV.transpose() * F_inv.transpose());
    SIGMAf_old=
      rho_f * nu_f * (F_inv_old * NV_old + NV_old.transpose() * F_inv_old.transpose());
  }
}

template <int DIM>
void DIV_proj<DIM>::Form(VectorIterator      b,
                         const FemFunction&  U,
                         const TestFunction& N) const {
  if (domain < 0)  // fluid
  {
    VECTOR phi;
    Multiplex::init_test<DIM>(phi, N);

    auto F_inv= F.inverse();

    // (rho_f*J*(F^(-T):grad(v), phi)
    b[0]+= rho_f * J * (F_inv.transpose().array() * NV.array()).sum() * N.m();

    // (J*sigma*F^-T, grad(phi))
    VECTOR X= J * (SIGMAf * F_inv.transpose()) * phi;
    for (int i= 0; i < DIM; ++i)
      b[i + 1]+= X(i, 0);

    X= -U[0].m() * J * F.inverse().transpose() * phi;
    for (int i= 0; i < DIM; ++i)
      b[i + 1]+= X(i, 0);

    // Old pressure
    X= (*OLD)[0].m() * J * F.inverse().transpose() * phi;
    for (int i= 0; i < DIM; ++i)
      b[i + 1]+= X(i, 0);

    // - (J_old*sigma_old*F_old^-T, grad(phi))
    X= J_old * (SIGMAf_old * F_inv.transpose()) * phi;
    for (int i= 0; i < DIM; ++i)
      b[i + 1]-= X(i, 0);
  }
}

template <int DIM>
void DIV_proj<DIM>::Matrix(EntryMatrix&        A,
                           const FemFunction&  U,
                           const TestFunction& M,
                           const TestFunction& N) const {
  if (domain < 0)  // fluid
  {
    VECTOR phi;
    Multiplex::init_test<DIM>(phi, N);
    VECTOR psi;
    Multiplex::init_test<DIM>(psi, M);
    //////////////// divergence
    //	b[0] += rho_f * J * (F.inverse().transpose().array() * NV.array()).sum() * N.m();
    //  wrt v
    for (int j= 0; j < DIM; ++j) {
      // wrt v
      A(0, j + 1)+= DIVERGENCE_V[j] * N.m();
/*
            // wrt u
            A(0, j + 1 + DIM) += rho_f * Jj[j] * divergence * N.m();

            A(0, j + 1 + DIM) += DIVERGENCE_U[j] * N.m();
*/        }

///////// tensor
// wrt V
for (int j= 0; j < DIM; ++j) {
  VECTOR X= TENSOR_dV[j] * phi;
  for (int i= 0; i < DIM; ++i)
    A(i + 1, j + 1)+= X(i, 0);
}
/*        // wrt U
        for (int j = 0; j < DIM; ++j)
        {
            VECTOR X = TENSOR_dU[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += X(i, 0);
        }
*/
// wrt P
VECTOR X= PRESSURE_P * phi;
for (int i= 0; i < DIM; ++i)
  A(i + 1, 0)+= X(i, 0);
/*        // wrt U
        for (int j = 0; j < DIM; ++j)
        {
            X = PRESSURE_U[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += X(i, 0);
        }

*/    }
}

template <int DIM>
void DIV_proj<DIM>::point_M(int j, const FemFunction& U, const TestFunction& M) const {
  if (domain < 0) {
    VECTOR psi;
    auto   F_inv= F.inverse();
    Multiplex::init_test<DIM>(psi, M);
    for (int j= 0; j < DIM; ++j) {
      // set to zero
      Jj[j] = (psi.transpose() * J * F_inv.block(0, j, DIM, 1))(0, 0);
      Fij[j]= -F_inv.block(0, j, DIM, 1) * psi.transpose() * F_inv;
    }
    divergence= (F_inv.transpose().array() * NV.array()).sum();

    for (int j= 0; j < DIM; ++j) {
      TENSOR_dV[j]= rho_f * nu_f * J
                    * (F_inv.block(0, j, DIM, 1) * psi.transpose()
                       + (F_inv.block(0, j, DIM, 1) * psi.transpose()).transpose())
                    * F_inv.transpose();
      TENSOR_dU[j]= rho_f * nu_f * J * (Fij[j] * NV + (Fij[j] * NV).transpose())
                      * F_inv.transpose()                  // p1
                    + J * SIGMAf * Fij[j].transpose()      // p2
                    + Jj[j] * SIGMAf * F_inv.transpose();  // p3

      DIVERGENCE_U[j]= rho_f * J * (Fij[j].transpose().array() * NV.array()).sum();
      DIVERGENCE_V[j]= rho_f * J * (psi.transpose() * F_inv.block(0, j, DIM, 1))(0, 0);

      PRESSURE_P= -M.m() * J * F_inv.transpose();

      PRESSURE_U[j]=
        -U[0].m() * Jj[j] * F_inv.transpose() - U[0].m() * J * Fij[j].transpose();
    }
  }
}

template <int DIM>
void DIV_proj<DIM>::MatrixBlock(EntryMatrix&       A,
                                const FemFunction& U,
                                const FemFunction& N) const {
  for (int j= 0; j < N.size(); ++j)  // trial
  {
#define M N[j]
    point_M(j, U, M);

    for (int i= 0; i < N.size(); ++i)  // test
    {
      A.SetDofIndex(i, j);
      Matrix(A, U, M, N[i]);
    }
#undef M
  }
}

////////////////////////////////////////////////// LPS

template <int DIM>
void DIV_proj<DIM>::lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const {
  double vel= 1.0;

  lps= lps0 / (vel / h + nu_f / h / h);
  // domain = -1;
  domain= chi(v);
}
template <int DIM>
void DIV_proj<DIM>::StabForm(VectorIterator      b,
                             const FemFunction&  U,
                             const FemFunction&  UP,
                             const TestFunction& N) const {
  if (domain < 0)  /// fludi
    for (int i= 0; i < DIM; ++i)
      b[0]+= lps * UP[0][i + 1] * N[i + 1];
}

template <int DIM>
void DIV_proj<DIM>::StabMatrix(EntryMatrix&        A,
                               const FemFunction&  U,
                               const TestFunction& Np,
                               const TestFunction& Mp) const {
  if (domain < 0)
    for (int i= 0; i < DIM; ++i)
      A(0, 0)+= lps * Mp[i + 1] * Np[i + 1];
}

template class DIV_proj<2>;
template class DIV_proj<3>;

}  // namespace Gascoigne

/*-----------------------------------------*/
