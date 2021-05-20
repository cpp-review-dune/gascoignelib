#include "boundarysolideuler.h"
#include "filescanner.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {

template<int DIM>
BoundarySolidEuler<DIM>::BoundarySolidEuler(const ParamFile* pf)
{
  DataFormatHandler DFH;
  DFH.insert("nu_f", &__nu_f, 0.0);
  DFH.insert("rho_f", &__rho_f);
  FileScanner FS(DFH, pf, "Equation");
  p_2 = 2.266 * 1.0e4;
  p_4 = 2.286 * 1.0e4;
  cout << "%%%%%%%%%%BoundarySolidEuler%%%%%%%%%%" << endl;
  cout << "Boundary 2/4 -- do-nothing with p=" << p_2
       << "g/cm/s^2=" << p_2 / 1333.22 << "mmHg and p=" << p_4
       << "g/cm/s^2=" << p_4 / 1333.22 << "mmHg" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
}

template<int DIM>
void
BoundarySolidEuler<DIM>::Form(VectorIterator b,
                              const FemFunction& U,
                              const TestFunction& N,
                              int col) const
{
  if (DIM == 3) {
    NV << (*VEL)[1].x(), (*VEL)[1].y(), (*VEL)[1].z(), (*VEL)[2].x(),
      (*VEL)[2].y(), (*VEL)[2].z(), (*VEL)[3].x(), (*VEL)[3].y(), (*VEL)[3].z();
  } else {
    NV << (*VEL)[1].x(), (*VEL)[1].y(), (*VEL)[2].x(), (*VEL)[2].y();
  }

  g = -__rho_f * __nu_f * NV.transpose() * __n;

  if (col == 4) {
    for (int i = 0; i < DIM; ++i) {
      b[i + 1] += g(i) * N.m() + __n[i] * p_2 * N.m();
      if (((*VEL)[1].m() * __n[0] + (*VEL)[2].m() * __n[1] +
           (*VEL)[3].m() * __n[2]) < 0)
        b[i + 1] -= 0.5 *
                    ((*VEL)[1].m() * __n[0] + (*VEL)[2].m() * __n[1] +
                     (*VEL)[3].m() * __n[2]) *
                    (*VEL)[1 + i].m() * N.m();
    }
  }

  if (col == 2) {
    for (int i = 0; i < DIM; ++i) {
      b[i + 1] += g(i) * N.m() + __n[i] * p_4 * N.m();
      if (((*VEL)[1].m() * __n[0] + (*VEL)[2].m() * __n[1] +
           (*VEL)[3].m() * __n[2]) < 0)
        b[i + 1] -= 0.5 *
                    ((*VEL)[1].m() * __n[0] + (*VEL)[2].m() * __n[1] +
                     (*VEL)[3].m() * __n[2]) *
                    (*VEL)[1 + i].m() * N.m();
    }
  }
}

template<int DIM>
void
BoundarySolidEuler<DIM>::Matrix(EntryMatrix& A,
                                const FemFunction& U,
                                const TestFunction& M,
                                const TestFunction& N,
                                int col) const
{}

template<int DIM>
void
BoundarySolidEuler<DIM>::pointboundary(double h,
                                       const FemFunction& U,
                                       const Vertex<DIM>& v,
                                       const Vertex<DIM>& n) const
{
  __n[0] = n[0];
  __n[1] = n[1];
  if (DIM == 3)
    __n[2] = n[2];
}

template class BoundarySolidEuler<2>;
template class BoundarySolidEuler<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
