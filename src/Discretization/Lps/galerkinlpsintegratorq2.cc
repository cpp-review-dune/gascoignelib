#include  "galerkinlpsintegratorq2.h"

namespace Gascoigne
{
/*-----------------------------------------*/

template<int DIM>
void GalerkinLpsIntegratorQ2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  GalerkinIntegrator<DIM>::Form(EQ,F,FEM,U,Q);
  Lps.                     Form(EQ,F,FEM,U,Q);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinLpsIntegratorQ2<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  GalerkinIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q);
  Lps.                     Matrix(EQ,E,FEM,U,Q);
}

/*-----------------------------------------------------------*/

template GalerkinLpsIntegratorQ2<2>;
template GalerkinLpsIntegratorQ2<3>;
}
