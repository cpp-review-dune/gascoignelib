#ifndef __IntegratorQ1Q2_h
#define __IntegratorQ1Q2_h

#include "basicintegrator.h"
#include "domainrighthandside.h"

namespace Gascoigne
{

/*---------------------------------------------------*/

template<int DIM>
class IntegratorQ1Q2 : public BasicIntegrator
{
 protected:

  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

 public:

  IntegratorQ1Q2<DIM>() : BasicIntegrator() {}
  ~IntegratorQ1Q2<DIM>() {}

  std::string GetName() const {return "IntegratorQ1Q2";}  

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, 
	   const FemInterface& FemL, const LocalNodeData& Q) const;
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
	    const FemInterface& FemL, const LocalVector& U, const LocalNodeData& Q) const;
  void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
	    const FemInterface& FemL, const LocalVector& U, const LocalNodeData& Q) const;
  void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int i, const LocalNodeData& Q) const;
};

}

/*---------------------------------------------------*/

#endif
