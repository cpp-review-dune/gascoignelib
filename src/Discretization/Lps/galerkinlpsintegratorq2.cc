#include  "galerkinlpsintegratorq2.h"

namespace Gascoigne
{
/*-----------------------------------------*/

template<int DIM>
void GalerkinLpsIntegratorQ2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegrator<DIM>::Form(EQ,F,FEM,U,Q,QC);
  Lps.                     Form(EQ,F,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinLpsIntegratorQ2<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q,QC);
  Lps.                     Matrix(EQ,E,FEM,U,Q,QC);
}

/*-----------------------------------------------------------*/

template class GalerkinLpsIntegratorQ2<2>;
template class GalerkinLpsIntegratorQ2<3>;
}
