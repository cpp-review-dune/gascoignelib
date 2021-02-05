#include "def_solid.h"
#include "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;
extern bool InterfaceResidual;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> Def_Solid<DIM>::Def_Solid(const ParamFile *pf) {}

//////////////////////////////////////////////////

template <int DIM> void Def_Solid<DIM>::point_cell(int material) const {
  if (material == 1)
    domain = 1;
  if (material == 2)
    domain = -1;
}

template <int DIM>
void Def_Solid<DIM>::point(double h, const FemFunction &U,
                           const Vertex<DIM> &v) const {
  __h = h;
  __v = v;
}

template <int DIM>
void Def_Solid<DIM>::Form(VectorIterator b, const FemFunction &U,
                          const TestFunction &N) const {
  if (domain > 0) // solid
  {
    for (int i = 0; i < DIM;
         i++) { //		(DEF-DEFOLD)-__THETA*VEL+(1-THETA)*VEL_OLD
      b[i] += (((*U_Vec)[DIM + 1 + i].m() - (*UOLD_Vec)[DIM + 1 + i].m()) -
               __DT * (__THETA * (*U_Vec)[i + 1].m() +
                       (1.0 - __THETA) * (*UOLD_Vec)[i + 1].m())) *
              N.m();
    }
  }
}

template <int DIM>
void Def_Solid<DIM>::Matrix(EntryMatrix &A, const FemFunction &U,
                            const TestFunction &M,
                            const TestFunction &N) const {
  if (domain > 0) {
    for (int i = 0; i < DIM; i++)
      A(i, i) += M.m() * N.m();
  }
}

template <int DIM>
void Def_Solid<DIM>::point_M(int j, const FemFunction &U,
                             const TestFunction &M) const {}

template <int DIM>
void Def_Solid<DIM>::MatrixBlock(EntryMatrix &A, const FemFunction &U,
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

////////////////////////////////////////////////// LPS

template <int DIM>
void Def_Solid<DIM>::lpspoint(double h, const FemFunction &U,
                              const Vertex<DIM> &v) const {}
template <int DIM>
void Def_Solid<DIM>::StabForm(VectorIterator b, const FemFunction &U,
                              const FemFunction &UP,
                              const TestFunction &N) const {}

template <int DIM>
void Def_Solid<DIM>::StabMatrix(EntryMatrix &A, const FemFunction &U,
                                const TestFunction &Np,
                                const TestFunction &Mp) const {}

template class Def_Solid<2>;
template class Def_Solid<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
