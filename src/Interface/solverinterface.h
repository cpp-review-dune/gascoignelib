#ifndef  __SolverInterface_h
#define  __SolverInterface_h

#include  "gascoigne.h"
#include  "matrixinterface.h"
#include  "iluinterface.h"

#include  "multigridmeshinterface.h"
#include  "functionalmanager.h"

#include  "problemdescriptorinterface.h"
#include  "meshinterface.h"
#include  "mginterpolatorinterface.h"
#include  "meshtransferinterface.h"
#include  "basicghostvector.h"
#include  "data.h"
#include  "paramfile.h"

//////////////////////////////////////////////
///
///@brief
/// Interface class for Solver

///  Some porporties
///  - lives on one level of the hierarchy
///  - stores the matrices
///  - holds memory for vectors
///  - provides e.g. nonlinear and linear residuals of the equations
///  - calls class MeshInterpretor
///
//////////////////////////////////////////////

/*---------------------------------------------------------*/

class SolverInterface
{
 protected:

 public:

  virtual ~SolverInterface(){}

  virtual std::string GetName() const=0;

  virtual void BasicInit(int level, const Gascoigne::ParamFile* paramfile, const MeshInterface* MP)=0;
  
  virtual void SetProblem(const ProblemDescriptorInterface& PDX)=0;
  virtual void SetState(const std::string& s)=0;

  virtual bool DirectSolver() const=0;
  
  virtual void NewMesh(int l, const MeshInterface* MP)=0;

  virtual void OutputSettings() const {assert(0);}
  virtual void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)=0;

  virtual void VisuGrid(const std::string& name, int i) const {assert(0);}

  virtual void AddNodeVector(const Gascoigne::GlobalVector* q) {assert(0);}
  virtual void AddCellVector(const Gascoigne::GlobalVector* q) {assert(0);}
  virtual void AddParameterVector(const Gascoigne::GlobalVector* q) {assert(0);}

  //
  /// vector - manamgement
  //

  virtual void ResizeVector(Gascoigne::GlobalVector* x, std::string type) const=0;
  virtual void RegisterVector(const BasicGhostVector& g)=0;
  virtual Gascoigne::GlobalVector& GetGV(BasicGhostVector& u) const=0;
  virtual const Gascoigne::GlobalVector& GetGV(const BasicGhostVector& u) const=0;

  //
  /// vector - hanging nodes
  //

  virtual void HNAverage   (const BasicGhostVector& x) const=0;
  virtual void HNZero      (const BasicGhostVector& x) const=0;
  virtual void HNDistribute(BasicGhostVector& x) const=0;

  virtual void HNAverage   (const Gascoigne::GlobalVector& x) const=0;
  virtual void HNZero      (const Gascoigne::GlobalVector& x) const=0;
  virtual bool HNZeroCheck(const Gascoigne::GlobalVector& x) const=0;
  virtual void HNDistribute(Gascoigne::GlobalVector& x) const=0;

  //
  /// vector - io
  //

  virtual void Visu(const std::string& name, const BasicGhostVector& u, int i) const=0;
  virtual void Write(const BasicGhostVector& u, const std::string& filename) const{assert(0);}
  virtual void Read(BasicGhostVector& u, const std::string& filename) const{assert(0);}

  //
  /// vector - interpolation
  //

  virtual void InterpolateSolution(BasicGhostVector& u, const Gascoigne::GlobalVector& uold) const=0;

  //
  /// vector - rhs (integration)
  //

  virtual void TimeRhs(BasicGhostVector& f, const BasicGhostVector& u) const {assert(0);}
  virtual void Rhs(BasicGhostVector& f, double d=1.) const=0;

  //
  /// vector - residual (integration)
  //

  virtual void Residual(BasicGhostVector& y, const BasicGhostVector& x, double d) const=0;

  //
  /// vector - boundary condition
  //

  virtual void SetBoundaryVector(Gascoigne::GlobalVector& f) const=0;
  virtual void SetBoundaryVector(BasicGhostVector& f) const=0;
  virtual void SetBoundaryVectorZero(BasicGhostVector& Gf) const=0;

  //
  /// vector - linear algebra
  //
  virtual double NewtonNorm(const BasicGhostVector& u) const=0;
  virtual void residualgmres(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const=0;
  virtual void MatrixResidual(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const=0;
  virtual void vmulteqgmres(BasicGhostVector& y, const BasicGhostVector& x) const=0;
  virtual void smooth_pre(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;
  virtual void smooth_exact(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;
  virtual void smooth_post(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const=0;

  //
  /// vector - additional
  //

  virtual void PressureFilterIntegrate(BasicGhostVector& x) const=0;

  //
  /// vector - matrix
  //

  virtual void AssembleMatrix(BasicGhostVector& u, double d=1.)=0;
  virtual void DirichletMatrix() const=0;
  virtual void MatrixZero() const=0;
  virtual void ComputeIlu(const BasicGhostVector& u) const=0;
  virtual void TransposeMatrix(const BasicGhostVector& u) const { assert(0);};

  //
  /// vector - "postprocessing"
  //

  virtual void ComputeError(const BasicGhostVector& u, Gascoigne::GlobalVector& err) const=0;
  virtual double ComputeFunctional(BasicGhostVector& f, const BasicGhostVector& u, const Functional* FP) const=0;
  virtual double EnergyEstimator(nvector<double>& eta, const BasicGhostVector& u, BasicGhostVector& f) const=0;

  //
  /// vector - initialize
  //

  virtual void BoundaryInit(BasicGhostVector& u) const=0;
  virtual void SolutionInit(BasicGhostVector& u) const=0;
};

#endif

