#include  "lpsintegrator.h"
#include  "patchintegrationformula.h"

namespace Gascoigne
{
/*----------------------------------------- */
/*   Q2 */
/*----------------------------------------- */

template<int DIM>
LpsIntegratorQ2<DIM>::LpsIntegratorQ2<DIM>() : LpsIntegrator<DIM>()
{
  if (DIM==2)
    _IF = new QuadGauss9;
  else
    _IF = new HexGauss27;
  assert(_IF);
}

/*----------------------------------------- */
/*   Q1 */
/*----------------------------------------- */

template<int DIM>
LpsIntegratorQ1<DIM>::LpsIntegratorQ1<DIM>() : LpsIntegrator<DIM>()
{
  if (DIM==2)
    {
      _IF = new PatchFormula2d<4,QuadGauss4>;
      CellWeight = 0.25;
    }
  else
    {
      _IF = new PatchFormula3d<8,HexGauss8>;
      CellWeight = 0.125;
    }
  assert(_IF);
}

/*-----------------------------------------*/

template<int DIM>
void LpsIntegratorQ1<DIM>::VectorReinit(LocalVector& F, int n, int ncomp) const
{
  F.ReInit(ncomp,n);
  F.zero();
}

/*-----------------------------------------*/

template<int DIM>
void LpsIntegratorQ1<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  VectorReinit(F,FEM.n(),EQ.ncomp());
  LpsIntegrator<DIM>::Form(EQ,F,FEM,U,Q);
}

/*-----------------------------------------*/

template<int DIM>
void LpsIntegratorQ1<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  ResetMatrix(E,FEM.n(),U.ncomp());
  LpsIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q);
}

/*-----------------------------------------------------------*/

template<int DIM>
void LpsIntegratorQ1<DIM>::ResetMatrix(EntryMatrix& E, int n, int ncomp) const
{
  E.SetDimensionDof(n,n);
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();
}

/*-----------------------------------------------------------*/

template LpsIntegratorQ1<2>;
template LpsIntegratorQ1<3>;
template LpsIntegratorQ2<2>;
template LpsIntegratorQ2<3>;
}
