#ifndef  __BechmarkFunctionals_h
#define  __BechmarkFunctionals_h

#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"

/*-----------------------------------------*/
/* Functionals for the 2d Navier-Stokes Benchmark*/
/*-----------------------------------------*/

namespace Gascoigne
{

class DragFunctional : public virtual Gascoigne::ResidualFunctional
{
  public:

  DragFunctional() : ResidualFunctional() 
    {
    __comps.push_back(1);
    __cols.insert(80);
    __scales.push_back(50);
    ExactValue() = 5.579535;       // fuer den runden
;

    __DD  = new Gascoigne::DirichletDataByColor(GetComps(),GetColors(),GetScales());
  }

  std::string GetName() const { return "DragFunctional";}
};

}

/*-----------------------------------------*/

#endif
