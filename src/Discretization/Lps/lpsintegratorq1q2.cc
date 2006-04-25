#include  "lpsintegrator.h"
#include  "patchintegrationformula.h"

namespace Gascoigne
{
/*----------------------------------------- */
/*   Q2 */
/*----------------------------------------- */

template<int DIM>
LpsIntegratorQ2<DIM>::LpsIntegratorQ2() : LpsIntegrator<DIM>()
{
  int femn;
  if (DIM==2)
    {
      LpsIntegrator<DIM>::_IF = new QuadGauss9;
      femn = 9;
    }
  else
    {
      LpsIntegrator<DIM>::_IF = new HexGauss27;
      femn = 27;
    }
  assert(LpsIntegrator<DIM>::_IF);

  LpsIntegrator<DIM>::MMM.resize(femn);
}

/*----------------------------------------- */
/*   Q1 */
/*----------------------------------------- */

template<int DIM>
LpsIntegratorQ1<DIM>::LpsIntegratorQ1() : LpsIntegrator<DIM>()
{
  int femn;
  if (DIM==2)
    {
      LpsIntegrator<DIM>::_IF = new PatchFormula2d<4,QuadGauss4>;
      LpsIntegrator<DIM>::CellWeight = 0.25;
      femn = 9;
    }
  else
    {
      LpsIntegrator<DIM>::_IF = new PatchFormula3d<8,HexGauss8>;
      LpsIntegrator<DIM>::CellWeight = 0.125;
      femn = 27;
    }
  assert(LpsIntegrator<DIM>::_IF);
  LpsIntegrator<DIM>::MMM.resize(femn);
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
void LpsIntegratorQ1<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  VectorReinit(F,FEM.n(),EQ.GetNcomp());
  LpsIntegrator<DIM>::Form(EQ,F,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template<int DIM>
void LpsIntegratorQ1<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  ResetMatrix(E,FEM.n(),U.ncomp());
  LpsIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q,QC);
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

template class LpsIntegratorQ1<2>;
template class LpsIntegratorQ1<3>;
template class LpsIntegratorQ2<2>;
template class LpsIntegratorQ2<3>;
}
