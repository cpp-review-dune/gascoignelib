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
  void BasicInit();

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, 
      const FemInterface& FemL, const LocalNodeData& Q, const LocalCellData& QC) const;
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
      const FemInterface& FemL, const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const;
  void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
      const FemInterface& FemL, const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const;
  void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
      int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const;
  void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
      const LocalVector& U, int ile, int col, LocalNodeData& Q, const LocalCellData& QC) const;
  void DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, const Vertex<DIM>& p, 
      const DiracRightHandSide& DRHS, int j, const LocalNodeData& Q, const LocalCellData& QC) const;

  double MassMatrix(EntryMatrix& E, const FemInterface& FemH, const FemInterface& FemL) const;
};

}

/*---------------------------------------------------*/

#endif
