#include "dwrlps2d.h"

namespace Gascoigne
{

/*-------------------------------------------------*/

void DwrLps2d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
//   nmatrix<double> T;
//   for(int iq=0; iq<GetPatchMesh()->npatches(); ++iq)
//     {
//       Transformation(T,iq);
//       GetFem()->ReInit(T);
      
//       GlobalToLocal(__U,u,iq);
//       GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__Q);
//       LocalToGlobal(f,__F,iq,d);
//     }
  Q1Lps2d::Form(f,u,EQ,d);
}

/*-------------------------------------------------*/

}
