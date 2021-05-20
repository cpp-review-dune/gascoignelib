#include "meshmotion.h"
#include "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;
extern bool InterfaceResidual;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {
template<int DIM>
MeshMotion<DIM>::MeshMotion(const ParamFile* pf)
{}

//////////////////////////////////////////////////

//#include "multiplex.xx"

template<int DIM>
void
MeshMotion<DIM>::point_cell(int material) const
{
  if (material == 1)
    domain = 1;
  if (material == 2)
    domain = -1;
}

template<int DIM>
void
MeshMotion<DIM>::point(double h,
                       const FemFunction& U,
                       const Vertex<DIM>& v) const
{
  // double dx = std::max(0.0,fabs(v.x()-0.4)-0.4);
  // double dy = std::max(0.0,fabs(v.y()-0.2)-0.01);
  // double dist = sqrt(dx*dx+dy*dy);
  // ext = 1.0 / (1.e-2+dist);
  ext = 1.0;
}

template<int DIM>
void
MeshMotion<DIM>::Form(VectorIterator b,
                      const FemFunction& U,
                      const TestFunction& N) const
{

  if (domain > 0)
    return;

  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      b[i] += ext * ((*U_Vec)[1 + DIM + i][j + 1] * N[j + 1] +
                     (1. - __THETA) / __THETA *
                       (*UOLD_Vec)[1 + DIM + i][j + 1] * N[j + 1]);
}

template<int DIM>
void
MeshMotion<DIM>::Matrix(EntryMatrix& A,
                        const FemFunction& U,
                        const TestFunction& M,
                        const TestFunction& N) const
{

  if (domain > 0)
    return;

  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      A(i, i) += ext * M[j + 1] * N[j + 1];
}

template<int DIM>
void
MeshMotion<DIM>::point_M(int j,
                         const FemFunction& U,
                         const TestFunction& M) const
{}

template<int DIM>
void
MeshMotion<DIM>::MatrixBlock(EntryMatrix& A,
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

////////////////////////////////////////////////// LPS

template<int DIM>
void
MeshMotion<DIM>::lpspoint(double h,
                          const FemFunction& U,
                          const Vertex<DIM>& v) const
{}
template<int DIM>
void
MeshMotion<DIM>::StabForm(VectorIterator b,
                          const FemFunction& U,
                          const FemFunction& UP,
                          const TestFunction& N) const
{}

template<int DIM>
void
MeshMotion<DIM>::StabMatrix(EntryMatrix& A,
                            const FemFunction& U,
                            const TestFunction& Np,
                            const TestFunction& Mp) const
{}

template class MeshMotion<2>;
template class MeshMotion<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
