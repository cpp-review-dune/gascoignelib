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
#include  "stopwatch.h"
#include  "hierarchicalmesh.h"
#include  "pressurefilter.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear StdSolver

///
///
//////////////////////////////////////////////

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


  // 3. Discretization

  DiscretizationInterface*    _ZP;

  // 4. Vektoren

  mutable GhostVectorAgent _NGVA;

  // 5. Anwendungsklassen

  const ProblemDescriptorInterface*      _PDX;

  // 6. Steuerparameter

  bool                _distribute;
  int                 _mylevel;

  mutable int         _ndirect;
  mutable bool        _directsolver;
  mutable std::string _discname;
  mutable std::string _matrixtype;

  SolverData          _Dat;
  mutable int         _PrimalSolve;
  const ParamFile*    _paramfile;

  // 5. sonstiges
  
  PressureFilter       _PF;
/*   double               omega_domain; */

  mutable StopWatch    _vm, _il, _so, _ca, _ci, _cs, _re;

  //
  //        Funktionen
  //

  // 0. Zugriff

  // 0.3 Matrizen

//   const MatrixInterface* GetMatrix() const {assert(_MAP); return _MAP;}
  MatrixInterface* GetMatrix() const { return _MAP;}
  MatrixInterface*& GetMatrixPointer() {assert(_MAP==NULL); return _MAP;}

//   const IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface*& GetIluPointer() {assert(_MIP==NULL); return _MIP;}
  
  virtual DiscretizationInterface*& GetDiscretizationPointer() {return _ZP;}

  // 1. Initialisierung 
	
	void SetDefaultValues(std::string discname, std::string matrixtype, int ndirect);

  virtual DiscretizationInterface* NewDiscretization(int dimension, const std::string& discname);

  virtual MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype); 
  virtual IluInterface* NewIlu(int ncomp, const std::string& matrixtype); 

  //
  /// new interface-function for indivisual size of vectors
  //

  virtual void Rhs(GlobalVector& f, double d=1.) const;

  virtual double ComputeFunctional(GlobalVector& f, const GlobalVector& u, const Functional* FP) const;
    
  virtual void smooth(int niter, GlobalVector& x, const GlobalVector& y, GlobalVector& h) const;
  virtual void PermutateIlu(const GlobalVector& u) const;
  virtual void modify_ilu(IluInterface& I,int ncomp) const;
  void Form(GlobalVector& y, const GlobalVector& x, double d) const;
  void AdjointForm(GlobalVector& y, const GlobalVector& x, double d) const;

  void MatrixResidual(GlobalVector& y, const GlobalVector& x, const GlobalVector& b) const;
  virtual void vmult(GlobalVector& y, const GlobalVector& x, double d) const;
  virtual void vmulteq(GlobalVector& y, const GlobalVector& x, double d) const;

  virtual DoubleVector IntegrateSolutionVector(const GlobalVector& u) const;
  virtual void _check_consistency(const Equation* EQ, const DiscretizationInterface* MP) const;
  virtual void DirichletMatrixOnlyRow() const;

  virtual void SetBoundaryVector(GlobalVector& f) const;
  virtual void SetBoundaryVectorZero(GlobalVector& f) const;
  virtual void SetBoundaryVectorStrong(GlobalVector& f, const BoundaryManager& BM, const DirichletData& DD) const;

 public:

  StdSolver();
  ~StdSolver();

  std::string GetName() const {return "StdSolver";}

  void BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP);
  void SetProblem(const ProblemDescriptorInterface& PDX);
  void SetDiscretization(DiscretizationInterface& DI, bool init=false);
  const ProblemDescriptorInterface* GetProblemDescriptor() const {assert(_PDX); return _PDX;}
  const ParamFile* GetParamfile() const { return _paramfile;}

  void NewMesh(int l, const MeshInterface* MP);

  const MeshInterface* GetMesh() const {return _MP;}

  // 0.2 Discretization

  const DiscretizationInterface* GetDiscretization() const {assert(_ZP); return _ZP;}
  DiscretizationInterface* GetDiscretization() {assert(_ZP); return _ZP;}

  void ReInitVector();
  void ReInitMatrix();

  virtual void SubtractMean(GlobalVector& gx) const;
  virtual void SubtractMeanAlgebraic(GlobalVector& gx) const;

  virtual double clock_vmult() const {return _vm.read();}
  virtual double clock_ilu  () const {return _il.read();}
  virtual double clock_solve() const {return _so.read();}
  virtual double clock_computematrix() const {return _ca.read();}
  virtual double clock_computeilu   () const {return _ci.read();}
  virtual double clock_computesolver() const {return _cs.read();}
  virtual double clock_residual     () const {return _re.read();}

  bool DirectSolver() const {return _directsolver;}

  void AddNodeVector(const std::string& name, const GlobalVector* q) {
    GetDiscretization()->AddNodeVector(name,q);
  }
  void AddCellVector(const std::string& name, const GlobalCellVector* q) {
    GetDiscretization()->AddCellVector(name,q);
  }
  void AddParameterVector(const std::string& name, const GlobalParameterVector* q) {
    GetDiscretization()->AddParameterVector(name,q);
  }
  void DeleteNodeVector(const std::string& name)  {
    GetDiscretization()->DeleteNodeVector(name);
  }
  void DeleteCellVector(const std::string& name) {
    GetDiscretization()->DeleteCellVector(name);
  }
  void DeleteParameterVector(const std::string& name) {
    GetDiscretization()->DeleteParameterVector(name);
  }

  void OutputSettings() const;
  virtual void Visu(const std::string& name, const GlobalVector& u, int i) const;

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void VisuGrid(const std::string& name, int i) const;
  
  //
  /// vector - manamgement
  //
    
  void RegisterMatrix();
  void ResizeVector(GlobalVector* x, std::string type) const;
  void RegisterVector(const BasicGhostVector& g);

        GlobalVector& GetGV(      BasicGhostVector& u) const { return _NGVA(u);}
  const GlobalVector& GetGV(const BasicGhostVector& u) const { return _NGVA(u);}

  //
  /// vector - hanging nodes
  //

  bool distribute() const { return _distribute; }
  void SetDistribute(bool dist) { _distribute = dist; }

  void HNAverage   (const BasicGhostVector& x) const;
  void HNZero      (const BasicGhostVector& x) const;
  void HNDistribute(BasicGhostVector& x) const;
  void HNAverage   (const GlobalVector& x) const;
  void HNZero      (const GlobalVector& x) const;
  bool HNZeroCheck(const GlobalVector& x) const;
  void HNDistribute(GlobalVector& x) const;
  void HNAverageData() const;
  void HNZeroData() const;

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

  void Form(BasicGhostVector& y, const BasicGhostVector& x, double d) const;

  //
  /// vector - boundary condition
  //
  void SetBoundaryVector(BasicGhostVector& f) const;
  void SetBoundaryVectorZero(BasicGhostVector& Gf) const;
  void SetBoundaryVectorStrong(BasicGhostVector& f, const BoundaryManager& BM, const DirichletData& DD) const;

  //
  /// vector - linear algebra
  //

  double NewtonNorm(const BasicGhostVector& u) const;
  void residualgmres(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const;
  void MatrixResidual(BasicGhostVector& y, const BasicGhostVector& x, const BasicGhostVector& b) const;
  void vmult  (BasicGhostVector& y, const BasicGhostVector& x, double d) const;
  void vmulteq(BasicGhostVector& y, const BasicGhostVector& x) const;
  void smooth_pre(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void smooth_exact(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void smooth_post(BasicGhostVector& y, const BasicGhostVector& x, BasicGhostVector& h) const;
  void Zero(BasicGhostVector& dst) const;

  //
  /// vector - additional
  //

  void SubtractMean(BasicGhostVector& x) const;
  void SubtractMeanAlgebraic(BasicGhostVector& x) const;

  //
  /// vector - matrix
  //

  void AssembleMatrix(const BasicGhostVector& u, double d);
  void DirichletMatrix() const;
  void MatrixZero() const;
  void ComputeIlu(const BasicGhostVector& u) const;
  void ComputeIlu() const;
  void AssembleDualMatrix(const BasicGhostVector& gu, double d);

  //
  /// vector - "postprocessing"
  //

  void ComputeError(const BasicGhostVector& u, GlobalVector& err) const;
  double ComputeFunctional(BasicGhostVector& f, const BasicGhostVector& u, const Functional* FP) const;

  virtual double ComputeBoundaryFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const BoundaryFunctional* FP) const;
  virtual double ComputeDomainFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const DomainFunctional* FP) const;
  virtual double ComputePointFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const PointFunctional* NFP) const;
  virtual double ComputeResidualFunctional(GlobalVector& f, const GlobalVector& u, GlobalVector& z, const ResidualFunctional* FP) const;

  //
  /// vector - initialize
  //

  void BoundaryInit(BasicGhostVector& u) const;
  void SolutionInit(BasicGhostVector& u) const;

  //
  /// HierarchicalMesh
  //

  virtual const HierarchicalMesh*& GetHierarchicalMeshPointer() { return _HM; }
  virtual const HierarchicalMesh*  GetHierarchicalMesh() const  { return _HM; }

  //
  /// for gmres
  //
  virtual void MemoryVector(BasicGhostVector& p);
  virtual void DeleteVector(BasicGhostVector* p) const;
  virtual double ScalarProduct(const BasicGhostVector& y, const BasicGhostVector& x) const;
  virtual void Equ(BasicGhostVector& dst, double s, const BasicGhostVector& src) const;
  virtual void Add(BasicGhostVector& dst, double s, const BasicGhostVector& src) const;
};
}

#endif



