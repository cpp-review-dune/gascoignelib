#include "galerkinglsintegrator.h"


/*-----------------------------------------*/

namespace Gascoigne
{
template<int DIM>
void GalerkinGlsIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegrator<DIM>::Form(EQ,F,FEM,U,Q,QC);
  Gls.                     Form(EQ,F,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q,QC);
  Gls.                     Matrix(EQ,E,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template class GalerkinGlsIntegrator<2>;
template class GalerkinGlsIntegrator<3>;
}
