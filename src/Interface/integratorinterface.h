#ifndef  __IntegratorInterface_h
#define  __IntegratorInterface_h


/////////////////////////////////////////////
///
///@brief
///  ... comments IntegratorInterface

///
///
/////////////////////////////////////////////

#include  "gascoigne.h"
#include  "equation.h"
#include  "feminterface.h"
#include  "entrymatrix.h"
#include  "domainrighthandside.h"
#include  "domainfunctional.h"
#include  "neumanndata.h"
#include  "robindata.h"

#include  "diracrighthandside.h"

namespace Gascoigne
{
class IntegratorInterface
{
public:

  IntegratorInterface() {}
  virtual ~IntegratorInterface() {}
  
  virtual std::string GetName() const=0;

  virtual void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const { assert(0);};
  virtual void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const { assert(0);};
  virtual void BoundaryForm(const RobinData& RD, LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile, int col, LocalNodeData& Q) const { assert(0); }
  virtual void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const { assert(0);};
  virtual void BoundaryMatrix (const RobinData& RD, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, int ile, int col, const LocalNodeData& Q) const { assert(0); }
  virtual double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const { assert(0); return 0;};

  virtual void RhsNeumann(const NeumannData& RHS, LocalVector& F, const FemInterface& FEM, int ile, int col, const LocalNodeData& Q) const {assert(0);}

  virtual double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const {assert(0); return 0;}

  virtual void ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, const LocalNodeData& Q) const  {assert(0);}

  virtual void RhsPoint(LocalVector& F, const FemInterface& FEM, const Vertex2d& v, int) const { assert(0);};
  virtual void RhsPoint(LocalVector& F, const FemInterface& FEM, const Vertex3d& v, int) const { assert(0);};
  virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex2d& p, const DiracRightHandSide& DRHS, int i, const LocalNodeData& Q) const { assert(0);};
  virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex3d& p, const DiracRightHandSide& DRHS, int i, const LocalNodeData& Q) const { assert(0);};
  virtual double ComputePointValue(const FemInterface& E, const Vertex2d& p, const LocalVector& U, int comp) const { assert(0);};
  virtual double ComputePointValue(const FemInterface& E, const Vertex3d& p, const LocalVector& U, int comp) const { assert(0);};
};
}

#endif
