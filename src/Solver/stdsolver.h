#ifndef  __StdSolver_h
#define  __StdSolver_h

#include  "gascoigne.h"
#include  "solverinterface.h"
#include  "gascoignemesh.h"
#include  "solverdata.h"
#include  "newghostvectoragent.h"
#include  "multigridmeshinterface.h"
#include  "pointfunctional.h"
#include  "residualfunctional.h"
#include  "meshinterpretorinterface.h"
#include  "stopwatch.h"

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear StdSolver

///
///
//////////////////////////////////////////////

using namespace Gascoigne;

/*-----------------------------------------*/

class StdSolver : public virtual SolverInterface
{
 private:

  //
  //   Daten
  //

  // 1. Gitter

  const MeshInterface*  _MP;

  // 2. Matrizen

  mutable MatrixInterface*  _MAP;
  mutable IluInterface*     _MIP;
  
 protected:

  // 3. MeshInterpretor

  MeshInterpretorInterface*    _ZP;

  // 4. Vektoren

  mutable NewGhostVectorAgent _NGVA;

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
  std::string         _paramfile;

  // 5. sonstiges
  
  nvector<double>      PF;
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

  const ProblemDescriptorInterface* GetProblemDescriptor() const {assert(_PDX); return _PDX;}

  // 0.3 Matrizen

//   const MatrixInterface* GetMatrix() const {assert(_MAP); return _MAP;}
  MatrixInterface* GetMatrix() const {assert(_MAP); return _MAP;}
  MatrixInterface*& GetMatrixPointer() {assert(_MAP==NULL); return _MAP;}

//   const IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface*& GetIluPointer() {assert(_MIP==NULL); return _MIP;}

  // 1. Initialisierung 

  virtual MeshInterpretorInterface* NewMeshInterpretor(int dimension, const std::string& discname);

  virtual MatrixInterface* NewMatrix(int ncomp, const string& matrixtype); 
  virtual IluInterface* NewIlu(int ncomp, const string& matrixtype); 

  virtual void MemoryMatrix();

  //
  /// new interface-function for indivisual size of vectors
  //

  std::string GetName() const {return "StdSolver";}

  void Rhs(GlobalVector& f, double d=1.) const;
  int RhsPoint(GlobalVector& f, const std::vector<Vertex2d>& p, int mycomp, const nvector<double>& s) const;

  double ComputeFunctional(GlobalVector& f, const GlobalVector& u, const Functional* FP) const;
  double ComputeBoundaryFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const BoundaryFunctional* FP) const;
  double ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const;
  double ComputePointFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const PointFunctional* FP) const;
  double ComputeResidualFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const ResidualFunctional* FP) const;
  
  void SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const Equation& EQ, const DirichletData& DD) const;
  virtual void smooth(int niter, GlobalVector& x, const GlobalVector& y, GlobalVector& h) const;
  void PressureFilter(GlobalVector& gx) const;
  void PressureFilterIntegrate(GlobalVector& gx) const;
  virtual void PermutateIlu(const GlobalVector& u) const;
  void modify_ilu(IluInterface& I,int ncomp) const;
  void ConstructPressureFilter();
  void Residual(GlobalVector& y, const GlobalVector& x, double d) const;

  void MatrixResidual(GlobalVector& y, const GlobalVector& x, const GlobalVector& b) const;
  void vmult(GlobalVector& y, const GlobalVector& x, double d) const;
  void vmulteq(GlobalVector& y, const GlobalVector& x, double d) const;

 public:

  StdSolver();
  ~StdSolver();

  void BasicInit(int level, const std::string& pfile, const MeshInterface* MP, const ProblemDescriptorInterface* PDX);

  void SetProblem(const ProblemDescriptorInterface& PDX);
  void NewMesh(int l, const MeshInterface* MP);

  virtual void MemoryVector();

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

  void AddNodeVector(const GlobalVector* q) {
    GetMeshInterpretor()->AddNodeVector(q);
  }
  void AddCellVector(const GlobalVector* q) {
    GetMeshInterpretor()->AddCellVector(q);
  }
  void AddParameterVector(const GlobalVector* q) {
    GetMeshInterpretor()->AddParameterVector(q);
  }
  void DeleteNodeVector(const GlobalVector* q)  {
    GetMeshInterpretor()->DeleteNodeVector(q);
  }
  void DeleteCellVector(const GlobalVector* q) {
    GetMeshInterpretor()->DeleteCellVector(q);
  }
  void DeleteParameterVector(const GlobalVector* q) {
    GetMeshInterpretor()->DeleteParameterVector(q);
  }

  void OutputSettings() const;

 void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void VisuGrid(const std::string& name, int i) const;

  //
  /// vector - manamgement
  //

  void ResizeVector(GlobalVector* x, string type) const;
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

  void Residual(BasicGhostVector& y, const BasicGhostVector& x, double d) const;

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

  //
  /// vector - additional
  //

  void PressureFilterIntegrate(BasicGhostVector& x) const;

  //
  /// vector - matrix
  //

  void AssembleMatrix(BasicGhostVector& u, double d);
  void DirichletMatrix() const;
  void MatrixZero() const;
  void ComputeIlu(const BasicGhostVector& u) const;

  //
  /// vector - "postprocessing"
  //

  void ComputeError(const BasicGhostVector& u, GlobalVector& err) const;
  double ComputeFunctional(BasicGhostVector& f, const BasicGhostVector& u, const Functional* FP) const;
  double EnergyEstimator(nvector<double>& eta, const BasicGhostVector& u, BasicGhostVector& f) const;

  //
  /// vector - initialize
  //

  void BoundaryInit(BasicGhostVector& u) const;
  void SolutionInit(BasicGhostVector& u) const;
};

#endif



