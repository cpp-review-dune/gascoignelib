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

class Q12d : public Q1
{
 protected:

  nmatrix<double> GetLocalInterpolationWeights() const;

  HNStructureInterface* NewHNStructure();

public:

  //
  ////  Con(De)structor 
  //
  
  Q12d();

  std::string GetName() const {return "Q12d";}
  
  void BasicInit(const Gascoigne::ParamFile* pf);

  void Interpolate(Gascoigne::GlobalVector& u, const InitialCondition& U) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(Gascoigne::GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const;


  void Jumps(EdgeInfoContainer<2>& EIC, const Gascoigne::GlobalVector& u) const;
  void JumpNorm(EdgeInfoContainer<2>& EIC, nvector<double>& eta) const;
};


#endif
