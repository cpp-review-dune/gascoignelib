#ifndef __DwrQ2Q4_h
#define __DwrQ2Q4_h

#include "solverinterface.h"

namespace Gascoigne
{

/**********************************************************/

  class DwrQ2Q4
  {
    private:
      SolverInterface                  &_S;
      const ProblemDescriptorInterface *_P;
      DiscretizationInterface          *_D;

    protected:
      DiscretizationInterface* CreateOtherDiscretization() const;

      double ScalarProduct(DoubleVector &eta, const GlobalVector &f, const GlobalVector &z) const;
      double ScalarProduct(DoubleVector &eta, const VectorInterface &gf, const VectorInterface &gz) const;
      double ScalarProductWithFluctuations(DoubleVector& eta, const VectorInterface &gf, const VectorInterface &gz) const;

      void PrimalResidualsHigher(VectorInterface &gf, const VectorInterface &gu);

      void DualResidualsHigher(VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI);

    public:
      DwrQ2Q4(SolverInterface &S);
      ~DwrQ2Q4() { }

      double Estimator(DoubleVector &eta, VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI);
  };

  /**********************************************************/

}

#endif
