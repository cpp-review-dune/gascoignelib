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
#include  "righthandsidedata.h"
#include  "domainfunctional.h"
#include  "neumanndata.h"

class IntegratorInterface
{
public:

  IntegratorInterface() {}
  virtual ~IntegratorInterface() {}
  
  virtual std::string GetName() const=0;

  virtual void Rhs(const RightHandSideData& RHS, Gascoigne::LocalVector& F, const FemInterface& FEM, const Gascoigne::LocalData& Q) const { assert(0);};
  virtual void Form(const Equation& EQ, Gascoigne::LocalVector& F, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const { assert(0);};
  virtual void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const { assert(0);};
  virtual void MassMatrix(EntryMatrix& E, const FemInterface& FEM) const { assert(0);};

  virtual void RhsNeumann(const NeumannData& RHS, Gascoigne::LocalVector& F, const FemInterface& FEM, int ile, int col, const Gascoigne::LocalData& Q) const {assert(0);}

  virtual double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const {assert(0);}

  virtual void ErrorsByExactSolution(Gascoigne::LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const Gascoigne::LocalVector& U, const Gascoigne::LocalData& Q) const  {assert(0);}

  virtual void RhsPoint(Gascoigne::LocalVector& F, const FemInterface& FEM, const Vertex2d& v, int) const { assert(0);};
  virtual void RhsPoint(Gascoigne::LocalVector& F, const FemInterface& FEM, const Vertex3d& v, int) const { assert(0);};
};


#endif
