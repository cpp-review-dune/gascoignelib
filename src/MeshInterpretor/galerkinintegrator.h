#ifndef  __GalerkinIntegrator_h
#define  __GalerkinIntegrator_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinIntegrator
////
////
/////////////////////////////////////////////

#include  "basicintegrator.h"
#include  "integrationformula.h"

template<int DIM>
class GalerkinIntegrator : public BasicIntegrator
{
private:

  IntegrationFormulaInterface*  IFF;
  IntegrationFormulaInterface*  IFE;
  IntegrationFormulaInterface*  IFB;

protected:

  IntegrationFormulaInterface*& FormFormulaPointer() { return IFF;}
  IntegrationFormulaInterface*& ErrorFormulaPointer() { return IFE;}
  IntegrationFormulaInterface*& BoundaryFormulaPointer() { return IFB;}

  const IntegrationFormulaInterface& FormFormula() const { return *IFF;}
  const IntegrationFormulaInterface& ErrorFormula() const { return *IFE;}
  const IntegrationFormulaInterface& BoundaryFormula() const { return *IFB;}

  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

public:

//
////  Con(De)structor 
//

  GalerkinIntegrator<DIM>();
  ~GalerkinIntegrator<DIM>() {}
  
  std::string GetName() const {return "Galerkin";}

  void Rhs(const RightHandSideData& RHS, Gascoigne::LocalVector& F, const FemInterface& FEM, const Gascoigne::LocalData& Q) const;
  void Form(const Equation& EQ, Gascoigne::LocalVector& F, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const;
  void MassMatrix(EntryMatrix& E, const FemInterface& FEM) const;

  void RhsPoint(Gascoigne::LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, int comp) const;
  double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const;

  void ErrorsByExactSolution(Gascoigne::LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const;

  void RhsNeumann(const NeumannData& RHS, Gascoigne::LocalVector& F, const FemInterface& FEM, int ile, int col, const Gascoigne::LocalData& Q) const;

};


#endif
