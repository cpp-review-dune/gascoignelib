#include "galerkinglsintegrator.h"

using namespace Gascoigne;

/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const
{
  GalerkinIntegrator<DIM>::Form(EQ,F,FEM,U,Q);
  Gls.                     Form(EQ,F,FEM,U,Q);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const
{
  GalerkinIntegrator<DIM>::Matrix(EQ,E,FEM,U,Q);
  Gls.                     Matrix(EQ,E,FEM,U,Q);
}

/*-----------------------------------------*/

template GalerkinGlsIntegrator<2>;
template GalerkinGlsIntegrator<3>;
