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

namespace Gascoigne
{
template<int DIM>
class GalerkinIntegrator : public BasicIntegrator
{
private:

  IntegrationFormulaInterface*  IFF;
  IntegrationFormulaInterface*  IFE;
  IntegrationFormulaInterface*  IFB;
  IntegrationFormulaInterface*  IFM;

protected:

  IntegrationFormulaInterface*& FormFormulaPointer() { return IFF;}
  IntegrationFormulaInterface*& ErrorFormulaPointer() { return IFE;}
  IntegrationFormulaInterface*& MassFormulaPointer() { return IFM;}
  IntegrationFormulaInterface*& BoundaryFormulaPointer() { return IFB;}

  const IntegrationFormulaInterface* FormFormula() const { assert(IFF); return IFF;}
  const IntegrationFormulaInterface* MassFormula() const { assert(IFM); return IFM;}
  const IntegrationFormulaInterface* ErrorFormula() const { assert(IFE); return IFE;}
  const IntegrationFormulaInterface* BoundaryFormula() const { assert(IFB); return IFB;}

  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

public:

//
////  Con(De)structor 
//

  GalerkinIntegrator<DIM>();
  ~GalerkinIntegrator<DIM>();
  
  std::string GetName() const {return "Galerkin";}
  void BasicInit();

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const;
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
  void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
  void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile, int col, LocalNodeData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
  void BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, int ile, int col, const LocalNodeData& Q) const;
  double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const;

  void RhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, int comp) const;
  void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalNodeData& Q) const;
  double ComputePointValue(const FemInterface& E, const Vertex<DIM>& p, const LocalVector& U, int comp) const;
  double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;

  void ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, const LocalNodeData& Q) const;

  void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, int ile, int col, const LocalNodeData& Q) const;
};
}


#endif
