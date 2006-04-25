#ifndef __IntegratorWithSecond_h
#define __IntegratorWithSecond_h

#include  "galerkinintegratorq2.h"
#include  "finiteelementwithsecond.h"
#include  "transformation2d.h"
#include  "baseq22dwithsecond.h"
#include  "transformation3d.h"
#include  "baseq23dwithsecond.h"
#include  "../Q1/finiteelement.xx"
#include  "finiteelementwithsecond.xx"
#include  "finiteelementwithsecond2d.xx"
#include  "finiteelementwithsecond3d.xx"

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
class IntegratorWithSecond : public GalerkinIntegratorQ2<DIM>
{
  protected:
 
  void point_hesse(const FemInterface& E, const Vertex<DIM>& v) const;

  void init_test_hesse(const FemInterface& E, TestFunction& N, double w, int i) const;

  void hesse(const FemInterface& E, FemFunction& UH, const LocalVector& u) const;

  void hesse(const FemInterface& E, FemData& QH, const LocalData& Q) const;

  public:
    
  std::string GetName() const {return "IntegratorWithSecond";}

  double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const; 

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, 
      const LocalData& Q, const LocalData& QC) const;

};
}

/**********************************************************/

#endif
