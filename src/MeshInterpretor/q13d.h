#ifndef  __Q1Simple3d_h
#define  __Q1Simple3d_h

#include  "q1.h"


#include  "edgeinfocontainer.h"

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Simple3d
////
////
/////////////////////////////////////////////

class Q13d : public Q1
{
 protected:

  nmatrix<double> GetLocalInterpolationWeights() const;

  HNStructureInterface* NewHNStructure();

public:

  //
  ////  Con(De)structor 
  //
  
  Q13d();

  std::string GetName() const {return "Q13d";}
  
  void BasicInit(const Gascoigne::ParamFile* pf);

  void Interpolate(Gascoigne::GlobalVector& u, const InitialCondition& U) const;
  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void StrongDirichletVector(Gascoigne::GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const;


  void Jumps(EdgeInfoContainer<3>&, const Gascoigne::GlobalVector&) const;
  void JumpNorm(EdgeInfoContainer<3>&, nvector<double>&) const;
};


#endif
