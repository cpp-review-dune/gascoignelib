#include "fluid_stat.h"
#include "filescanner.h"

using namespace std;

/*--------------------------------------------------------------------*/

namespace Gascoigne {
template<int DIM>
Fluid_Stat<DIM>::Fluid_Stat(const ParamFile* pf)
{
  DataFormatHandler DFH;

  DFH.insert("rho_f", &rho_f, 0.0);
  DFH.insert("nu_f", &nu_f, 0.0);
  DFH.insert("lps", &lps0, 0.0);
  FileScanner FS(DFH, pf, "Equation");
  assert(rho_f > 0);
  assert(nu_f > 0);
  assert(lps0 > 0);

  cout << "%%%%%%%%%% Problem Fluid stat %%%%%%%%%%" << endl
       << "  rho_f / nu_f: " << rho_f << " / " << nu_f << endl
       << "   lps: " << lps0 << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
}

//////////////////////////////////////////////////

#include "multiplex_fluid_stat.xx"

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::point_cell(int material) const
{
  if (material == 1)
    domain = 1;
  if (material == 2)
    domain = -1;
}
/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::point(double h,
                       const FemFunction& U,
                       const Vertex<DIM>& v) const
{
  __h = h;
  __v = v;

  if (domain < 0) {
    multiplex_fluid_stat_init_NV<DIM>(NV, U);
    multiplex_fluid_stat_init_V<DIM>(V, U);

    // TENSOR
    SIGMAf = rho_f * nu_f * (NV + NV.transpose());
  }
}

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::Form(VectorIterator b,
                      const FemFunction& U,
                      const TestFunction& N) const
{
  VECTOR phi;
  multiplex_fluid_stat_init_test<DIM>(phi, N);

  if (domain < 0) // fluid
  {
    // divergence
    b[0] += rho_f * NV.trace() * N.m();

    // tensor
    VECTOR X = SIGMAf * phi;
    for (int i = 0; i < DIM; ++i)
      b[i + 1] += X(i, 0);

    // convection
    X = rho_f * NV * V * N.m();
    for (int i = 0; i < DIM; ++i)
      b[i + 1] += X(i, 0);

    // pressure
    X = -U[0].m() * MATRIX::Identity() * phi;
    for (int i = 0; i < DIM; ++i)
      b[i + 1] += X(i, 0);
  }
}

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::Matrix(EntryMatrix& A,
                        const FemFunction& U,
                        const TestFunction& M,
                        const TestFunction& N) const
{

  VECTOR phi;
  multiplex_fluid_stat_init_test<DIM>(phi, N);
  VECTOR psi;
  multiplex_fluid_stat_init_test<DIM>(psi, M);

  if (domain < 0) // fluid
  {

    //////////////// divergence

    for (int j = 0; j < DIM; ++j) {
      // wrt v
      A(0, j + 1) += DIVERGENCE_V[j] * N.m();
    }

    ///////// tensor
    // wrt V
    for (int j = 0; j < DIM; ++j) {
      VECTOR X = TENSOR_dV[j] * phi;
      for (int i = 0; i < DIM; ++i)
        A(i + 1, j + 1) += X(i, 0);
    }

    // //////////////// // convection

    // wrt v
    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        A(i + 1, j + 1) += CONV_dV1(i, j) * M.m() * N.m();

    for (int j = 0; j < DIM; ++j)
      A(j + 1, j + 1) += CONV_dV2[j] * N.m();

    /////////////// pressure
    VECTOR X = PRESSURE_P * phi;
    for (int i = 0; i < DIM; ++i)
      A(i + 1, 0) += X(i, 0);
  }
}

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::point_M(int j,
                         const FemFunction& U,
                         const TestFunction& M) const
{
  VECTOR psi;
  multiplex_fluid_stat_init_test<DIM>(psi, M);

  if (domain < 0) {

    CONV_dV1 = rho_f * NV;

    for (int jj = 0; jj < DIM; ++jj) {
      TENSOR_dV[jj] =
        rho_f * nu_f *
        (MATRIX::Identity().block(0, jj, DIM, 1) * psi.transpose() +
         (MATRIX::Identity().block(0, jj, DIM, 1) * psi.transpose())
           .transpose());

      DIVERGENCE_V[jj] = rho_f * psi[jj];

      CONV_dV2[jj] = rho_f * (psi.transpose() * V)(0, 0);

      PRESSURE_P = -M.m() * MATRIX::Identity();
    }
  }
}

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::MatrixBlock(EntryMatrix& A,
                             const FemFunction& U,
                             const FemFunction& N) const
{
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

/*--------------------------------------------------------------------*/
////////////////////////////////////////////////// LPS

template<int DIM>
void
Fluid_Stat<DIM>::lpspoint(double h,
                          const FemFunction& U,
                          const Vertex<DIM>& v) const
{

  // double vel =10.0;
  double vel =
    1.0 * sqrt(U[1].m() * U[1].m() + U[2].m() * U[2].m() + U[3].m() * U[3].m());

  lps = lps0 / (vel / h + nu_f / h / h);
}
/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::StabForm(VectorIterator b,
                          const FemFunction& U,
                          const FemFunction& UP,
                          const TestFunction& N) const
{
  if (domain < 0) /// fludi
    for (int i = 0; i < DIM; ++i)
      b[0] += lps * UP[0][i + 1] * N[i + 1];
}

/*--------------------------------------------------------------------*/
template<int DIM>
void
Fluid_Stat<DIM>::StabMatrix(EntryMatrix& A,
                            const FemFunction& U,
                            const TestFunction& Np,
                            const TestFunction& Mp) const
{
  if (domain < 0)
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += lps * Mp[i + 1] * Np[i + 1];
}

template class Fluid_Stat<2>;
template class Fluid_Stat<3>;

} // namespace Gascoigne

/*--------------------------------------------------------------------*/
