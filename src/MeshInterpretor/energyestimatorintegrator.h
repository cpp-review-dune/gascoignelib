#ifndef __EnergyEstimatorIntegrator_h
#define __EnergyEstimatorIntegrator_h

#include "basicintegrator.h"
#include "integrationformula.h"

/**********************************************************/

template<int DIM>
class EnergyEstimatorIntegrator : public BasicIntegrator
{
  protected:
    const IntegrationFormulaInterface* IF;
    fixarray<2*DIM-2,Vertex<DIM-1> >   xi;

    double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}
    const IntegrationFormulaInterface& GetFormula() const {return *IF;}

  public:
    EnergyEstimatorIntegrator<DIM>() : BasicIntegrator() { }
    ~EnergyEstimatorIntegrator<DIM>();

    void BasicInit();
    std::string GetName() const {return "EnergyEstimator";}

    void   Jumps(Gascoigne::LocalVector& F, const FemInterface& FEM, const Gascoigne::LocalVector& U, int ile) const;
    double JumpNorm(const FemInterface& FEM, fixarray<2*DIM-2,double> jumps, int ile) const;
    double Residual(const Gascoigne::LocalVector& U, const FemInterface& FEM, const Equation& EQ, const RightHandSideData& RHS, const Gascoigne::LocalData& Q) const;
    double ResidualZeroRhs(const Gascoigne::LocalVector& U, const FemInterface& FEM, const Equation& EQ, const Gascoigne::LocalData& Q) const;
};

/**********************************************************/

#endif
