#ifndef  __StdSolver_h
#define  __StdSolver_h

#include  "gascoigne.h"
#include  "solverinterface.h"
#include  "gascoignemesh.h"
#include  "solverdata.h"
#include  "ghostvectoragent.h"
#include  "multigridmeshinterface.h"
#include  "pointfunctional.h"
#include  "residualfunctional.h"
#include  "meshinterpretorinterface.h"
#include  "stopwatch.h"
#include  "hierarchicalmesh.h"
#include  "pressurefilter.h"

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear StdSolver

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

namespace Gascoigne
{
class StdSolver : public virtual SolverInterface
{
 private:

  //
  //   Daten
  //

  // 1. Gitter

  const MeshInterface*    _MP;
  const HierarchicalMesh* _HM;

  // 2. Matrizen

  mutable MatrixInterface*  _MAP;
  mutable IluInterface*     _MIP;
  
 protected:


  // 3. MeshInterpretor

  MeshInterpretorInterface*    _ZP;

  // 4. Vektoren

  mutable GhostVectorAgent _NGVA;

  // 5. Anwendungsklassen

  const ProblemDescriptorInterface*      _PDX;

  // 6. Steuerparameter

  int                 _mylevel;

  mutable int         _ndirect;
  mutable bool        _directsolver;
  mutable std::string _discname;
  mutable std::string _matrixtype;

  SolverData          Dat;
  mutable int         _PrimalSolve;
  const ParamFile*    _paramfile;

  // 5. sonstiges
  
  PressureFilter       PF;
  double               omega_domain;

  mutable StopWatch    _vm, _il, _so, _ca, _ci, _cs, _re;

  //
  //        Funktionen
  //

  // 0. Zugriff

  // 0.1 Gitter

  const MeshInterface* GetMesh() const {return _MP;}

  // 0.2 MeshInterpretor

  MeshInterpretorInterface*& GetMeshInterpretorPointer() {assert(_ZP==NULL); return _ZP;}
  const MeshInterpretorInterface* GetMeshInterpretor() const {assert(_ZP); return _ZP;}
  MeshInterpretorInterface* GetMeshInterpretor() {assert(_ZP); return _ZP;}

  // 0.3 Matrizen

//   const MatrixInterface* GetMatrix() const {assert(_MAP); return _MAP;}
  MatrixInterface* GetMatrix() const {assert(_MAP); return _MAP;}
  MatrixInterface*& GetMatrixPointer() {assert(_MAP==NULL); return _MAP;}

//   const IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface*& GetIluPointer() {assert(_MIP==NULL); return _MIP;}

  // 1. Initialisierung 

  virtual MeshInterpretorInterface* NewMeshInterpretor(int dimension, const std::string& discname);

  virtual MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype); 
  virtual IluInterface* NewIlu(int ncomp, const std::string& matrixtype); 

  //
  /// new interface-function for indivisual size of vectors
  //

  std::string GetName() const {return "StdSolver";}

  void Rhs(GlobalVector& f, double d=1.) const;
  int RhsPoint(GlobalVector& f, const PointFunctional* FP) const;

  double ComputeFunctional(GlobalVector& f, const GlobalVector& u, const Functional* FP) const;
  double ComputeBoundaryFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const BoundaryFunctional* FP) const;
  double ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const;
  double ComputePointFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const PointFunctional* FP) const;
  
  void SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const DirichletData& DD) const;
  virtual void smooth(int niter, GlobalVector& x, const GlobalVector& y, GlobalVector& h) const;
  void SubtractMean(GlobalVector& gx) const;
  void SubtractMeanAlgebraic(GlobalVector& gx) const;
  virtual void PermutateIlu(const GlobalVector& u) const;
  void modify_ilu(IluInterface& I,int ncomp) const;
  void Form(GlobalVector& y, const GlobalVector& x, double d) const;

  void MatrixResidual(GlobalVector& y, const GlobalVector& x, const GlobalVector& b) const;
  void vmult(GlobalVector& y, const GlobalVector& x, double d) const;
  void vmulteq(GlobalVector& y, const GlobalVector& x, double d) const;

  DoubleVector IntegrateSolutionVector(const GlobalVector& u) const;
  virtual void _check_consistency(const Equation* EQ, const MeshInterpretorInterface* MP) const;

 public:

  StdSolver();
  ~StdSolver();

  void BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP);
  void SetProblem(const ProblemDescriptorInterface& PDX);
  const ProblemDescriptorInterface* GetProblemDescriptor() const {assert(_PDX); return _PDX;}

  void NewMesh(int l, const MeshInterface* MP);

  void ReInitVector();
  void ReInitMatrix();

  double clock_vmult() const {return _vm.read();}
  double clock_ilu  () const {return _il.read();}
  double clock_solve() const {return _so.read();}
  double clock_computematrix() const {return _ca.read();}
  double clock_computeilu   () const {return _ci.read();}
  double clock_computesolver() const {return _cs.read();}
  double clock_residual     () const {return _re.read();}

  void SetState(const std::string& s) {
    if     (s=="State")   _PrimalSolve = 1;
    else if(s=="Adjoint") _PrimalSolve = 0;
    else if(s=="Tangent") _PrimalSolve = 2;
    else abort();
  }
  
  bool DirectSolver() const {return _directsolver;}

  void AddNodeVector(const std::string& name, const GlobalVector* q) {
    HNAverage(*q);
    GetMeshInterpretor()->AddNodeVector(name,q);
  }
  void AddCellVector(const std::string& name, const GlobalCellVector* q) {
    GetMeshInterpretor()->AddCellVector(name,q);
  }
  void AddParameterVector(const std::string& name, const GlobalParameterVector* q) {
    GetMeshInterpretor()->AddParameterVector(name,q);
  }
  void DeleteNodeVector(const std::string& name)  {
    GetMeshInterpretor()->DeleteNodeVector(name);
  }
  void DeleteCellVector(const std::string& name) {
    GetMeshInterpretor()->DeleteCellVector(name);
  }
  void DeleteParameterVector(const std::string& name) {
    GetMeshInterpretor()->DeleteParameterVector(name);
  }

  void OutputSettings() const;

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void VisuGrid(const std::string& name, int i) const;
  
  //
  /// vector - manamgement
  //
    
  void RegisterMatrix();
  void ResizeVector(GlobalVector* x, std::string type) const;
  void RegisterVector(const BasicGhostVector& g) {_NGVA.Register(g,this);}
  GlobalVector& GetGV(BasicGhostVector& u) const {
    return _NGVA(u);
  }
  const GlobalVector& GetGV(const BasicGhostVector& u) const {
    return _NGVA(u);
  }

  //
  /// vector - hanging nodes
  //

  void HNAverage   (const BasicGhostVector& x) const;
  void HNZero      (const BasicGhostVector& x) const;
  void HNDistribute(BasicGhostVector& x) const;
  void HNAverage   (const GlobalVector& x) const;
  void HNZero      (const GlobalVector& x) const;
  bool HNZeroCheck(const GlobalVector& x) const;
  void HNDistribute(GlobalVector& x) const;

  //
  /// vector - io
  //

  void Visu(const std::string& name, const BasicGhostVector& u, int i) const;
  void Visu(const std::string& name, const GlobalVector& u, int i) const;
  void Write(const BasicGhostVector& u, const std::string& filename) const;
  void Read(BasicGhostVector& u, const std::string& filename) const;

  //
  /// vector - interpolation
  //

  void InterpolateSolution(BasicGhostVector& u, const GlobalVector& uold) const;

  //
  /// vector - rhs (integration)
  //

  void Rhs(BasicGhostVector& f, double d=1.) const;

  //
  /// vector - residual (integration)
  //

  void Form(BasicGhostVector& y, const BasicGhostVector& x, double d) const;

  //
  /// vector - boundary condition
  //

  void SetBoundaryVector(GlobalVector& f) const;
  void SetBoundaryVector(BasicGhostVector& f) const;
  void SetBoundaryVectorZero(BasicGhostVector& Gf) const;
  void SetBoundaryVectorZero(GlobalVector& f) const;

  //
  /// vector - linear algebra
  //

  double NewtonNorm(const BasicGhostVector& u) const;
  void residualgmres(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const;
  void MatrixResidual(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const;
  void vmulteqgmres(BasicGhostVector& y, const BasicGhostVector& x) const;
  void smooth_pre(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void smooth_exact(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void smooth_post(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void Zero(BasicGhostVector& dst) const;

  //
  /// vector - additional
  //

  void SubtractMean(BasicGhostVector& x) const;

  //
  /// vector - matrix
  //

  void AssembleMatrix(const BasicGhostVector& u, double d);
  void DirichletMatrix() const;
  void MatrixZero() const;
  void ComputeIlu(const BasicGhostVector& u) const;
  void ComputeIlu() const;

  //
  /// vector - "postprocessing"
  //

  void ComputeError(const BasicGhostVector& u, GlobalVector& err) const;
  double ComputeFunctional(BasicGhostVector& f, const BasicGhostVector& u, const Functional* FP) const;
  double EnergyEstimator(DoubleVector& eta, const BasicGhostVector& u, BasicGhostVector& f) const;

  //
  /// vector - initialize
  //

  void BoundaryInit(BasicGhostVector& u) const;
  void SolutionInit(BasicGhostVector& u) const;

  //
  /// HierarchicalMesh
  //

  const HierarchicalMesh*& GetHierarchicalMeshPointer() { return _HM; }
  const HierarchicalMesh*  GetHierarchicalMesh() const  { return _HM; }


  virtual double ComputeResidualFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const ResidualFunctional* FP) const;
};
}

#endif



