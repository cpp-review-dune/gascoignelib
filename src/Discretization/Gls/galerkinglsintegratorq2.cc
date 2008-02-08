#include  "galerkinglsintegratorq2.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegratorQ2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegratorQ2<DIM>::Form(EQ,F,FEM,U,Q,QC);
  Gls.                       Form(EQ,F,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegratorQ2<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegratorQ2<DIM>::Matrix(EQ,E,FEM,U,Q,QC);
  Gls.                       Matrix(EQ,E,FEM,U,Q,QC);
}

/*-----------------------------------------------------------*/

#if __GNUC__ == 3
template GalerkinGlsIntegratorQ2<2>;
template GalerkinGlsIntegratorQ2<3>;
#endif
#if __GNUC__ == 4
template class GalerkinGlsIntegratorQ2<2>;
template class GalerkinGlsIntegratorQ2<3>;
#endif
}
