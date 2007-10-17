#ifndef  __Q1Simple_h
#define  __Q1Simple_h

#include  "q1.h"
#include  "edgeinfocontainer.h"
#include  "energyestimatorintegrator.h"
#include  "curve.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Simple
////
////
/////////////////////////////////////////////

class Q12d : public Q1
{
 protected:

  HNStructureInterface* NewHNStructure();

  void EEJumps(EdgeInfoContainer<2>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEJumpNorm(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const;
  void EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const EnergyEstimatorIntegrator<2>& EEI) const;
  int GetCellNumber(const Vertex2d& p0, Vertex2d& p,int c0=0) const;
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
  void InterpolateDirac(GlobalVector& u, const GlobalVector& uold)const;

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const;

  void EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype, double d_visc) const;

  nmatrix<double> GetLocalInterpolationWeights() const;

  void RhsCurve(GlobalVector &F, const Curve &C,int comp,int N) const;
  
  Vertex2d randpunkt(const Curve& C, double t0,double& t1,int& cc) const;
};
}

#endif
