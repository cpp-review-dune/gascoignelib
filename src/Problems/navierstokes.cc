#include "navierstokes.h"
#include "filescanner.h"

namespace Gascoigne
{
  ////////////////////////////////////////////////// DATA

  void NavierStokesData::BasicInit(const ParamFile *pf)
  {
    DataFormatHandler DFH;
    DFH.insert("visc", &visc, 1.);
    DFH.insert("alpha", &alpha0, 0.);
    DFH.insert("dt", &dt, 0.);
    FileScanner FS(DFH, pf, "Equation");
    assert(visc>0);
  }

  ////////////////////////////////////////////////// NavierStokes

  template <int DIM>
  void NavierStokes<DIM>::Form(VectorIterator b,
                               const FemFunction &U,
                               const TestFunction &N) const
  {
    for (int i = 0; i < DIM; ++i)
    {
      b[0] += U[i + 1][i + 1] * N.m();
      for (int j = 0; j < DIM; ++j)
      {

        b[i + 1] += data.visc * U[i + 1][j + 1] * N[j + 1];
        b[i + 1] += U[j + 1].m() * U[i + 1][j + 1] * N.m();
      }

      b[i + 1] -= U[0].m() * N[i + 1];
    }
  }

  template <int DIM>
  void NavierStokes<DIM>::Matrix(EntryMatrix &A,
                                 const FemFunction &U,
                                 const TestFunction &M,
                                 const TestFunction &N) const
  {
    for (int i = 0; i < DIM; ++i)
    {
      A(0, i + 1) += M[i + 1] * N.m();
      for (int j = 0; j < DIM; ++j)
      {

        A(i + 1, i + 1) += data.visc * M[j + 1] * N[j + 1];

        A(i + 1, j + 1) += M.m() * U[i + 1][j + 1] * N.m();
        A(i + 1, i + 1) += U[j + 1].m() * M[j + 1] * N.m();
      }

      A(i + 1, 0) -= M.m() * N[i + 1];
    }
  }

  ////////////////////////////////////////////////// NavierStokesLps
  
  template<int DIM>
  void NavierStokesLps<DIM>::lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const
  {
    _h = h;
    _alpha = NavierStokes<DIM>::data.alpha0 * _h * _h/ NavierStokes<DIM>::data.visc;
  }
  
  template<int DIM>
  void NavierStokesLps<DIM>::StabForm(VectorIterator b,
				      const FemFunction &U,
				      const FemFunction &UP,
				      const TestFunction &NP) const
  {
    for (int i = 0; i < DIM; ++i)
      b[0] += _alpha * UP[0][i + 1] * NP[i + 1];
  }
  
  template<int DIM>
  void NavierStokesLps<DIM>::StabMatrix(EntryMatrix &A,
					const FemFunction &U,
					const TestFunction &Np,
					const TestFunction &Mp) const
  {
    for (int i = 0; i < DIM; ++i)
      A(0, 0) += _alpha * Mp[i + 1] * Np[i + 1];
  }
  
  
  template class NavierStokes<2>;
  template class NavierStokes<3>;

  template class NavierStokesLps<2>;
  template class NavierStokesLps<3>;

  
} // namespace Gascoigne
