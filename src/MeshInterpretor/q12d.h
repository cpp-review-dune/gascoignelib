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

class Q12d : public Q1
{
 protected:

  nmatrix<double> GetLocalInterpolationWeights() const;

  HNStructureInterface* NewHNStructure();

  void EEJumps(EdgeInfoContainer<2>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEJumpNorm(EdgeInfoContainer<2>& EIC, nvector<double>& eta, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEResidual(nvector<double>& eta, const GlobalVector& u, const Equation& EQ, const RightHandSideData& RHS, const EnergyEstimatorIntegrator<2>& EEI) const;
  void EEResidualZeroRhs(nvector<double>& eta, const GlobalVector& u, const Equation& EQ, const EnergyEstimatorIntegrator<2>& EEI) const;

public:

  //
  ////  Con(De)structor 
  //
  
  Q12d();

  std::string GetName() const {return "Q12d";}
  
  void BasicInit(const Gascoigne::ParamFile* pf);

  void Interpolate(GlobalVector& u, const InitialCondition& U) const;
  void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const;

  void EnergyEstimator(EdgeInfoContainer<2>& EIC, nvector<double>& eta, const GlobalVector& u, const Equation& EQ, const RightHandSideData& RHS) const;
  void EnergyEstimatorZeroRhs(EdgeInfoContainer<2>& EIC, nvector<double>& eta, const GlobalVector& u, const Equation& EQ) const;
};


#endif
