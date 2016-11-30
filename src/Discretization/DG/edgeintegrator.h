#ifndef __edgeintegrator__h
#define __edgeintegrator__h

#include "edgeequation.h"
#include "gascoigne.h"
#include "dgedge.h"


namespace Gascoigne
{

  class EdgeIntegratorBase
  {
  protected:
    //    IntegrationFormulaInterface*  EF;

  public:
    EdgeIntegratorBase() {}
    virtual ~EdgeIntegratorBase() {}
    
  };

  template<int DIM>
    class EdgeIntegrator : virtual public EdgeIntegratorBase
    {

    public:
      EdgeIntegrator<DIM>(){}

      //
      double Volume2MeshSize(const double V) const
      {
	return pow(V,1.0/static_cast<double>(DIM));
      }
      

      void EdgeForm(const EdgeEquation& EQ, const DGEdge& edge,
		    LocalVector& F1, LocalVector& F2,
		    const FemInterface& FEM1, const FemInterface& FEM2,
		    const LocalVector& U1, const LocalVector& U2,
		    const LocalData& Q, const LocalData& QC) const;

    };
  
}

#endif
