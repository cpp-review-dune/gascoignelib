#ifndef  __MeshInterpretorInterface_h
#define  __MeshInterpretorInterface_h


/////////////////////////////////////////////
///
///@brief
///  ... comments MeshInterpretorInterface

///
///
/////////////////////////////////////////////

#include  <string>
#include  "meshinterface.h"
#include  "gascoigne.h"
#include  "equation.h"
#include  "matrixinterface.h"
#include  "dirichletdata.h"
#include  "righthandsidedata.h"
#include  "neumanndata.h"
#include  "initialcondition.h"
#include  "exactsolution.h"
#include  "boundaryfunctional.h"
#include  "domainfunctional.h"
#include  "mginterpolatorinterface.h"
#include  "meshtransferinterface.h"
#include  "globaldata.h"
#include  "paramfile.h"

class MeshInterpretorInterface
{
private:

  mutable Gascoigne::GlobalData __q;

protected:

  const Gascoigne::GlobalData& GetGlobalData() const {return __q;}

public:

  MeshInterpretorInterface() {}
  virtual ~MeshInterpretorInterface() {}

  //
  //// Functions called from the Solver
  //

  virtual std::string GetName() const=0;

  virtual void AddNodeVector(const Gascoigne::GlobalVector* q) const {__q.AddNodeVector(q);}
  virtual void DeleteNodeVector(const Gascoigne::GlobalVector* q) const {__q.DeleteNodeVector(q);}
  virtual void AddCellVector(const Gascoigne::GlobalVector* q) const {__q.AddCellVector(q);}
  virtual void DeleteCellVector(const Gascoigne::GlobalVector* q) const {__q.DeleteCellVector(q);}
  virtual void AddParameterVector(const Gascoigne::GlobalVector* q) const {__q.AddParameterVector(q);}
  virtual void DeleteParameterVector(const Gascoigne::GlobalVector* q) const {__q.DeleteParameterVector(q);}

  virtual void BasicInit(const Gascoigne::ParamFile* pf)=0;
  virtual void ReInit   (const MeshInterface* M) { assert(0);};

  virtual int n() const=0;
  virtual int n_withouthanging()const {return n();}

  virtual void Structure(SparseStructureInterface* S) const { assert(0);};
  virtual void Form(Gascoigne::GlobalVector& f, const Gascoigne::GlobalVector& u, const Equation& EQ, double d) const { assert(0);};
  virtual void Matrix(MatrixInterface& A, const Gascoigne::GlobalVector& u, const Equation& EQ, double) const { assert(0);};
  virtual void MassMatrix(MatrixInterface& M) const {assert(0);}
  virtual void Rhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const { assert(0);};
  virtual void DiracRhs(Gascoigne::GlobalVector& f, const RightHandSideData& RHS, double s) const { assert(0);}
  virtual int RhsPoint(Gascoigne::GlobalVector& f, const Functional* F) const { assert(0);}
  virtual void RhsNeumann(Gascoigne::GlobalVector& f, const Equation& EQ, const Gascoigne::IntSet& Colors,  const NeumannData& NRHS, double s) const { assert(0);}


  virtual void HNAverage   (Gascoigne::GlobalVector& x) const {}
  virtual void HNDistribute(Gascoigne::GlobalVector& x) const {}
  virtual void HNZero      (Gascoigne::GlobalVector& x) const {}
  virtual bool HNZeroCheck (const Gascoigne::GlobalVector& x) const { return 0;}
  virtual void Interpolate(Gascoigne::GlobalVector& u, const InitialCondition& RHS) const { assert(0);};
  virtual void InterpolateSolution(Gascoigne::GlobalVector& u, const Gascoigne::GlobalVector& uold)const { assert(0);};
  virtual void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const { assert(0);};
  virtual void StrongDirichletVector(Gascoigne::GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const { assert(0);};
  virtual void StrongDirichletVectorZero(Gascoigne::GlobalVector& u, int col, const std::vector<int>& comp) const { assert(0);};


  virtual double PressureFilter(nvector<double>&) const {assert(0);}

  // Functionals
  virtual void ComputeError(const Gascoigne::GlobalVector& u, Gascoigne::LocalVector& err, const ExactSolution* ES) const { assert(0);};
  virtual double ComputeBoundaryFunctional(const Gascoigne::GlobalVector& u, const BoundaryFunctional& BF) const {assert(0);}
  virtual double ComputeDomainFunctional(const Gascoigne::GlobalVector& u, const DomainFunctional& F) const {assert(0);}

  virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT) { assert(0);};
};


#endif
