#ifndef  __Q1Simple3d_h
#define  __Q1Simple3d_h

#include  "q1.h"
#include  "edgeinfocontainer.h"
#include  "energyestimatorintegrator.h"

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Simple3d
////
////
/////////////////////////////////////////////

namespace Gascoigne
{
class Q13d : public Q1
{
 protected:

  nmatrix<double> GetLocalInterpolationWeights() const;

  HNStructureInterface* NewHNStructure();
  
  void EEJumps(EdgeInfoContainer<3>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const;
  void EEJumpNorm(EdgeInfoContainer<3>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const;
  void EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS, const EnergyEstimatorIntegrator<3>& EEI) const;
  void EEResidualZeroRhs(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const EnergyEstimatorIntegrator<3>& EEI) const;
  int GetCellNumber(const Vertex3d& p0, Vertex3d& p) const;
  void VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const;

public:

  //
  ////  Con(De)structor 
  //
  
  Q13d();

  std::string GetName() const {return "Q13d";}
  
  void BasicInit(const ParamFile* pf);

  void Interpolate(GlobalVector& u, const InitialCondition& U) const;
  void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const;

  void EnergyEstimator(EdgeInfoContainer<3>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS) const;
  void EnergyEstimatorZeroRhs(EdgeInfoContainer<3>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ) const;
};
}

#endif
