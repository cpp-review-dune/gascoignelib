#ifndef __EnergyEstimatorIntegrator_h
#define __EnergyEstimatorIntegrator_h

#include "basicintegrator.h"
#include "integrationformula.h"

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
class EnergyEstimatorIntegrator : public BasicIntegrator
{
  protected:
    const IntegrationFormulaInterface* IF;
    fixarray<2*DIM-2,Vertex<DIM-1> >   _xi;
    std::string _s_energytype;
    double      _d_visc;


    double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}
    const IntegrationFormulaInterface& GetFormula() const {return *IF;}

  public:

    EnergyEstimatorIntegrator<DIM>();
    EnergyEstimatorIntegrator<DIM>(const std::string & s_energytype,double d_visc);
    ~EnergyEstimatorIntegrator<DIM>();

    void BasicInit();
    std::string GetName() const {return "EnergyEstimator";}

    void   Jumps(LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile) const;
    double JumpNorm(const FemInterface& FEM, fixarray<2*DIM-2,double> jumps, int ile) const;
    double Residual(const LocalVector& U, const FemInterface& FEM, const Equation& EQ, const DomainRightHandSide* RHS, const LocalNodeData& Q) const;
};
}

/**********************************************************/

#endif
