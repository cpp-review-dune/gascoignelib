
/*----------------------------   dgintegrator.h ---------------------------*/
/*      $Id:$                 */
#ifndef __dgintegrator_H
#define __dgintegrator_H
/*----------------------------   dgintegrator.h ---------------------------*/

#include "dgequation.h"
#include "galerkinintegrator.h"

namespace Gascoigne
{


  template <int DIM>
  class DGIntegrator : public GalerkinIntegrator<DIM>
  {
  public:

    // Constructor / Init
    DGIntegrator();

    // Integration on edges
    void EdgeForm(bool internaledge,
		  const DGEquation &EQ,
                  LocalVector &F1,
                  LocalVector &F2,
                  const FemInterface &FEMASTER,
                  const FemInterface &FESLAVE,
		  const int& masterli,
		  const int& slaveli,
                  const LocalVector &U1,
                  const LocalVector &U2,
                  const LocalData &Q,
                  const LocalData &QC) const;
    void EdgeMatrix(bool internaledge,
		    const DGEquation &EQ,
		    EntryMatrix &E11,
		    EntryMatrix &E12,
		    EntryMatrix &E21,
		    EntryMatrix &E22,
		    const FemInterface &FEMASTER,
		    const FemInterface &FESLAVE,
		    const int& masterli,
		    const int& slaveli,
		    const LocalVector &U1,
		    const LocalVector &U2,
		    const LocalData &Q,
		    const LocalData &QC) const;
  };
}


/*----------------------------   dgintegrator.h ---------------------------*/
/* end of #ifndef __dgintegrator_H */
#endif
/*----------------------------   dgintegrator.h ---------------------------*/
