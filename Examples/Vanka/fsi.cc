#include "fsi.h"
#include "filescanner.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> EQ<DIM>::EQ(const ParamFile *pf) {
  DataFormatHandler DFH;
  DFH.insert("mu", &mu, 0.0);
  DFH.insert("lambda", &lambda, 0.0);
  FileScanner FS(DFH, pf, "Equation");
  assert(lambda >= 0);
  assert(mu > 0);
}

template <int DIM>
void EQ<DIM>::Form(VectorIterator b, const FemFunction &U,
                   const TestFunction &N) const {
  b[1] += N.m();
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j) {
      b[i] += mu * (U[i][j + 1] + U[j][i + 1]) * N[j + 1];
      b[i] += lambda * U[j][j + 1] * N[i + 1];
    }
}

template <int DIM>
void EQ<DIM>::Matrix(EntryMatrix &A, const FemFunction &U,
                     const TestFunction &M, const TestFunction &N) const {
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j) {
      A(i, i) += mu * M[j + 1] * N[j + 1];
      A(i, j) += mu * M[i + 1] * N[j + 1];
      A(i, j) += lambda * M[j + 1] * N[i + 1];
    }
}

template class EQ<2>;
template class EQ<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
