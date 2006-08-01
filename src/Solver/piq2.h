#ifndef __PiQ2_h
#define __PiQ2_h

#include "discretizationinterface.h"
#include "gascoignemesh.h"
#include "solverinterface.h"

namespace Gascoigne
{

/**********************************************************/

  class PiQ2
  {
    protected:
      const DiscretizationInterface *_DI;
      const GascoigneMesh           *_MP;
      nvector<DoubleVector>          _q2weight;

    public:
      PiQ2();
      ~PiQ2() { }

      void Init(const SolverInterface* SI);
      void vmult(GlobalVector& y, const GlobalVector& x) const;
  };

/**********************************************************/

}

#endif
