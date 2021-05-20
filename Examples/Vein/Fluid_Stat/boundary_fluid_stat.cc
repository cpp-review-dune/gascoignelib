#include "boundary_fluid_stat.h"
#include "filescanner.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {

////////////////////////////////////////////////// BOUNDARY
template<int DIM>
Boundary_Fluid_Stat<DIM>::Boundary_Fluid_Stat(const ParamFile* pf)
{
  DataFormatHandler DFH;
  DFH.insert("nu_f", &__nu_f, 0.0);
  DFH.insert("rho_f", &__rho_f);
  FileScanner FS(DFH, pf, "Equation");
  p_2 = 2.266 * 1.0e4;
  p_4 = 2.286 * 1.0e4;
  cout << "%%%%%%%%%%Fluid_Stat%%%%%%%%%%" << endl;
  cout << "Boundary 2/4 -- do-nothing with p=" << p_2
       << "g/cm/s^2=" << p_2 / 1333.22 << "mmHg and p=" << p_4
       << "g/cm/s^2=" << p_4 / 1333.22 << "mmHg" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
}

template<int DIM>
void
Boundary_Fluid_Stat<DIM>::Form(VectorIterator b,
                               const FemFunction& U,
                               const TestFunction& N,
                               int col) const
{
  //______________________________________________________________________________
  if (DIM == 3) {
    NV << U[1].x(), U[1].y(), U[1].z(), U[2].x(), U[2].y(), U[2].z(), U[3].x(),
      U[3].y(), U[3].z();
  } else {
    NV << U[1].x(), U[1].y(), U[2].x(), U[2].y();
  }

  g = -__rho_f * __nu_f * NV.transpose() * __n;

  if (col == 4) {
    for (int i = 0; i < DIM; ++i) {
      b[i + 1] += g(i) * N.m() + __n[i] * p_2 * N.m();
      if ((U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) < 0)
        b[i + 1] -=
          0.5 * (U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) *
          U[1 + i].m() * N.m();
    }
  }

  if (col == 2) {
    for (int i = 0; i < DIM; ++i) {
      b[i + 1] += g(i) * N.m() + __n[i] * p_4 * N.m();
      if ((U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) < 0)
        b[i + 1] -=
          0.5 * (U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) *
          U[1 + i].m() * N.m();
    }
  }
}

template<int DIM>
void
Boundary_Fluid_Stat<DIM>::Matrix(EntryMatrix& A,
                                 const FemFunction& U,
                                 const TestFunction& M,
                                 const TestFunction& N,
                                 int col) const
{

  //_______________________________________________________________
  if (DIM == 3) {
    NV << U[1].x(), U[1].y(), U[1].z(), U[2].x(), U[2].y(), U[2].z(), U[3].x(),
      U[3].y(), U[3].z();
  } else {
    NV << U[1].x(), U[1].y(), U[2].x(), U[2].y();
  }

  //________________________________________________________________

  for (int j = 0; j < DIM; ++j) {
    NPHI = MATRIX::Zero();
    PHI = VECTOR::Zero();
    for (int i = 0; i < DIM; ++i) {
      NPHI(j, i) = M[i + 1];
      PHI(j) = M.m();
    }
    //_________________________________________
    g = -__rho_f * __nu_f * NPHI.transpose() * __n;

    //___________________________________________
    for (int i = 0; i < DIM; ++i) {
      A(1 + i, j + 1) += g(i) * N.m();

      if ((U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) < 0) {
        A(1 + i, j + 1) -=
          0.5 * (U[1].m() * __n[0] + U[2].m() * __n[1] + U[3].m() * __n[2]) *
          PHI[i] * N.m();
        A(1 + i, j + 1) -=
          0.5 * (PHI[0] * __n[0] + PHI[1] * __n[1] + PHI[2] * __n[2]) *
          U[1 + i].m() * N.m();
      }
    }
  }
}

template<int DIM>
void
Boundary_Fluid_Stat<DIM>::pointboundary(double h,
                                        const FemFunction& U,
                                        const Vertex<DIM>& v,
                                        const Vertex<DIM>& n) const
{

  __n[0] = n[0];
  __n[1] = n[1];
  if (DIM == 3)
    __n[2] = n[2];
}

template class Boundary_Fluid_Stat<2>;
template class Boundary_Fluid_Stat<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
