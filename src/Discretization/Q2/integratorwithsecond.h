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

template<int DIM>
class IntegratorWithSecond : public Gascoigne::GalerkinIntegratorQ2<DIM>
{
  protected:
 
  void point_hesse(const FemInterface& E, const Gascoigne::Vertex<DIM>& v) const;

  void init_test_hesse(const Gascoigne::FemInterface& E, Gascoigne::TestFunction& N, double w, int i) const;

  void hesse(const Gascoigne::FemInterface& E, Gascoigne::FemFunction& UH, const Gascoigne::LocalVector& u) const;

  void hesse(const Gascoigne::FemInterface& E, Gascoigne::FemData& QH, const LocalNodeData& Q) const;

  public:
    
  std::string GetName() const {return "IntegratorWithSecond";}

  double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const; 

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const;

};

/**********************************************************/

#endif
