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

using namespace Gascoigne;


class MeshInterpretorInterface
{
private:

  mutable GlobalData __q;

protected:

  const GlobalData& GetGlobalData() const {return __q;}

public:

  MeshInterpretorInterface() {}
  virtual ~MeshInterpretorInterface() {}

  //
  //// Functions called from the Solver
  //

  virtual std::string GetName() const=0;

  virtual void AddNodeVector(const GlobalVector* q) const {__q.AddNodeVector(q);}
  virtual void DeleteNodeVector(const GlobalVector* q) const {__q.DeleteNodeVector(q);}
  virtual void AddCellVector(const GlobalVector* q) const {__q.AddCellVector(q);}
  virtual void DeleteCellVector(const GlobalVector* q) const {__q.DeleteCellVector(q);}
  virtual void AddParameterVector(const GlobalVector* q) const {__q.AddParameterVector(q);}
  virtual void DeleteParameterVector(const GlobalVector* q) const {__q.DeleteParameterVector(q);}

  virtual void BasicInit(const std::string& paramfile)=0;
  virtual void ReInit   (const MeshInterface* M) { assert(0);};

  virtual int n() const=0;
  virtual int n_withouthanging()const {return n();}

  virtual void Structure(SparseStructureInterface* S) const { assert(0);};
  virtual void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const { assert(0);};
  virtual void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double) const { assert(0);};
  virtual void MassMatrix(MatrixInterface& M) const {assert(0);}
  virtual void Rhs(GlobalVector& f, const RightHandSideData& RHS, double s) const { assert(0);};
  virtual void DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const { assert(0);}
  virtual int RhsPoint(GlobalVector& f, const std::vector<Vertex2d>& p0, int comp, const nvector<double>& d) const { assert(0);}
  virtual void RhsNeumann(GlobalVector& f, const Equation& EQ, const IntSet& Colors,  const NeumannData& NRHS, double s) const { assert(0);}


  virtual void HNAverage   (GlobalVector& x) const {}
  virtual void HNDistribute(GlobalVector& x) const {}
  virtual void HNZero      (GlobalVector& x) const {}
  virtual bool HNZeroCheck (const GlobalVector& x) const { return 0;}
  virtual void Interpolate(GlobalVector& u, const InitialCondition& RHS) const { assert(0);};
  virtual void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const { assert(0);};
  virtual void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const { assert(0);};
  virtual void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const { assert(0);};
  virtual void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const { assert(0);};


  virtual double PressureFilter(nvector<double>&) const {assert(0);}

  // Functionals
  virtual void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const { assert(0);};
  virtual double ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const {assert(0);}
  virtual double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const {assert(0);}

  virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT) { assert(0);};
};


#endif
