#ifndef  __Q1Simple_h
#define  __Q1Simple_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Simple
////
////
/////////////////////////////////////////////

#include  "q1.h"
#include  "edgeinfocontainer.h"
#include  "energyestimatorintegrator.h"

namespace Gascoigne
{
class Q12d : public Q1
{
 protected:

  nmatrix<double> GetLocalInterpolationWeights() const;

  HNStructureInterface* NewHNStructure();

  void EEJumps(EdgeInfoContainer<2>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEJumpNorm(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS, const EnergyEstimatorIntegrator<2>& EEI) const;
  void EEResidualZeroRhs(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const EnergyEstimatorIntegrator<2>& EEI) const;
  int GetCellNumber(const Vertex2d& p0, Vertex2d& p) const;
  void VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const;

public:

  //
  ////  Con(De)structor 
  //
  
  Q12d();

  std::string GetName() const {return "Q12d";}
  
  void BasicInit(const ParamFile* pf);

  void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const;
  void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const;

  void EnergyEstimator(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS) const;
  void EnergyEstimatorZeroRhs(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ) const;
};
}

#endif