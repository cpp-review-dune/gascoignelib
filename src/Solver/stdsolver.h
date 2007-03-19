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
#include  "domainfunction.h"
#include  "facediscretization.h"

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

  MatrixInterface*  _MAP;
  IluInterface*     _MIP;
  
 protected:


  // 3. Discretization

  DiscretizationInterface*    _ZP;
  FaceDiscretization*         _FZP;

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
  mutable std::string _facediscname;
  mutable std::string _matrixtype;

  SolverData          _Dat;
  mutable int         _PrimalSolve;
  const ParamFile*    _paramfile;

  bool _useUMFPACK;

  // 5. sonstiges
  
  PressureFilter       _PF;
/*   double               omega_domain; */

  mutable StopWatch    _vm, _il, _so, _ca, _ci, _cs, _re;

  //
  //        Funktionen
  //

  // 0. Zugriff

  const MeshInterface*& GetMeshPointer() {return _MP;}
  
  // 0.3 Matrizen

  MatrixInterface* GetMatrix() const { return _MAP;}
  MatrixInterface*& GetMatrixPointer() {return _MAP;}

  IluInterface* GetIlu() const {assert(_MIP); return _MIP;}
  IluInterface*& GetIluPointer() { return _MIP;}
  
  virtual DiscretizationInterface*& GetDiscretizationPointer()     {return _ZP;}
  virtual FaceDiscretization*&      GetFaceDiscretizationPointer() {return _FZP;}

  // 1. Initialisierung 
	
	void SetDefaultValues(std::string discname, std::string matrixtype, int ndirect);

  virtual DiscretizationInterface* NewDiscretization    (int dimension, const std::string& discname);
  virtual FaceDiscretization*      NewFaceDiscretization(int dimension, const std::string& facediscname);

  virtual MatrixInterface* NewMatrix(int ncomp, const std::string& matrixtype); 
  virtual IluInterface* NewIlu(int ncomp, const std::string& matrixtype); 

  //
  /// new interface-function for individual size of vectors
  //

  virtual void smooth(int niter, VectorInterface& x, const VectorInterface& y, VectorInterface& h) const;
  virtual void PermutateIlu(const VectorInterface& gu) const;
  virtual void modify_ilu(IluInterface& I,int ncomp) const;

  virtual DoubleVector IntegrateSolutionVector(const VectorInterface& u) const;
  virtual void _check_consistency(const Equation* EQ, const DiscretizationInterface* DI) const;
  virtual void DirichletMatrixOnlyRow() const;

 public:

  StdSolver();
  ~StdSolver();

  std::string GetName() const {return "StdSolver";}

  void BasicInit(int level, const ParamFile* paramfile, const int dimension);
  void SetProblem(const ProblemDescriptorInterface& PDX);
  void SetDiscretization(DiscretizationInterface& DI, bool init=false);
  const ProblemDescriptorInterface* GetProblemDescriptor() const {assert(_PDX); return _PDX;}
  const ParamFile* GetParamfile() const { return _paramfile;}

  void NewMesh(int l, const MeshInterface* MP);

  const MeshInterface* GetMesh() const {return _MP;}

  // 0.2 Discretization

  const DiscretizationInterface* GetDiscretization() const {assert(_ZP); return _ZP;}
  DiscretizationInterface* GetDiscretization() {assert(_ZP); return _ZP;}

  const FaceDiscretization* GetFaceDiscretization() const { return _FZP; }
  FaceDiscretization*       GetFaceDiscretization()       { return _FZP; }

  void ReInitMatrix();

  virtual double clock_vmult() const {return _vm.read();}
  virtual double clock_ilu  () const {return _il.read();}
  virtual double clock_solve() const {return _so.read();}
  virtual double clock_computematrix() const {return _ca.read();}
  virtual double clock_computeilu   () const {return _ci.read();}
  virtual double clock_computesolver() const {return _cs.read();}
  virtual double clock_residual     () const {return _re.read();}

  bool DirectSolver() const {return _directsolver;}

  void AddNodeVector(const std::string& name, const VectorInterface& q) {
    assert(q.GetType()=="node");
    GetDiscretization()->AddNodeVector(name,&GetGV(q));
  }
  void AddCellVector(const std::string& name, const VectorInterface& q) {
    assert(q.GetType()=="cell");
    GetDiscretization()->AddCellVector(name,&GetGV(q));
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
  virtual void PointVisu(const std::string& name, const GlobalVector& u, int i) const;
  virtual void CellVisu(const std::string& name, const GlobalVector& u, int i) const;

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);
  void VisuGrid(const std::string& name, int i) const;
  
  //
  /// vector - manamgement
  //
    
  void RegisterMatrix();
  void RegisterVector(const VectorInterface& g);
  void ReInitVector(VectorInterface& dst);
  void ReInitVector(VectorInterface& dst, int comp);

        GlobalVector& GetGV(      VectorInterface& u) const { return _NGVA(u);}
  const GlobalVector& GetGV(const VectorInterface& u) const { return _NGVA(u);}

  //
  /// vector - hanging nodes
  //

  bool GetDistribute() const { return _distribute; }
  void SetDistribute(bool dist) { _distribute = dist; }

  void HNAverage   (const VectorInterface& x) const;
  void HNZero      (const VectorInterface& x) const;
  void HNDistribute(VectorInterface& x) const;
  void HNAverageData() const;
  void HNZeroData() const;

  //
  /// vector - io
  //

  void Visu(const std::string& name, const VectorInterface& u, int i) const;
  void Write(const VectorInterface& u, const std::string& filename) const;
  void Read(VectorInterface& u, const std::string& filename) const;

  //
  /// vector - interpolation
  //

  void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const;

  //
  /// vector - rhs (integration)
  //

  void Rhs(VectorInterface& f, double d=1.) const;

  //
  /// vector - residual (integration)
  //

  void Form(VectorInterface& y, const VectorInterface& x, double d) const;
  void AdjointForm(VectorInterface& y, const VectorInterface& x, double d) const;

  //
  /// vector - boundary condition
  //
  void SetBoundaryVector(VectorInterface& f) const;
  void SetBoundaryVectorZero(VectorInterface& Gf) const;
  void SetBoundaryVectorStrong(VectorInterface& f, const BoundaryManager& BM, const DirichletData& DD, double d=1.) const;

  //
  /// vector - linear algebra
  //

  double NewtonNorm(const VectorInterface& u) const;
  void residualgmres(VectorInterface& y, const VectorInterface& x, const VectorInterface& b) const;
  void MatrixResidual(VectorInterface& y, const VectorInterface& x, const VectorInterface& b) const;
  void vmult  (VectorInterface& y, const VectorInterface& x, double d) const;
  void vmulteq(VectorInterface& y, const VectorInterface& x, double d) const;
  void smooth_pre(VectorInterface& y, const VectorInterface& x, VectorInterface& h) const;
  void smooth_exact(VectorInterface& y, const VectorInterface& x, VectorInterface& h) const;
  void smooth_post(VectorInterface& y, const VectorInterface& x, VectorInterface& h) const;
  void Zero(VectorInterface& dst) const;

  //
  /// vector - additional
  //

  void SubtractMean(VectorInterface& x) const;
  void SubtractMeanAlgebraic(VectorInterface& x) const;

  //
  /// vector - matrix
  //

  void AssembleMatrix(const VectorInterface& u, double d);
  void DirichletMatrix() const;
  void MatrixZero() const;
  void ComputeIlu(const VectorInterface& u) const;
  void ComputeIlu() const;
  void AssembleDualMatrix(const VectorInterface& gu, double d);

  //
  /// vector - "postprocessing"
  //

  void ComputeError(const VectorInterface& u, GlobalVector& err) const;
  void AssembleError(GlobalVector& eta, const VectorInterface& u, GlobalVector& err) const;
  double ComputeFunctional(VectorInterface& f, const VectorInterface& u, const Functional* FP) const;

  virtual double ComputeBoundaryFunctional(VectorInterface& f, const VectorInterface& u, VectorInterface& z, const BoundaryFunctional* FP) const;
  virtual double ComputeDomainFunctional(VectorInterface& f, const VectorInterface& u, VectorInterface& z, const DomainFunctional* FP) const;
  virtual double ComputePointFunctional(VectorInterface& f, const VectorInterface& u, VectorInterface& z, const PointFunctional* NFP) const;
  virtual double ComputeResidualFunctional(VectorInterface& f, const VectorInterface& u, VectorInterface& z, const ResidualFunctional* FP) const;
  virtual void EvaluateCellRightHandSide(VectorInterface& f, const DomainRightHandSide& CF, double d = 1.) const;
  virtual void InterpolateDomainFunction(VectorInterface& f, const DomainFunction& DF) const;

  //
  /// vector - initialize
  //

  void BoundaryInit(VectorInterface& u) const;
  void SolutionInit(VectorInterface& u) const;

  //
  /// HierarchicalMesh
  //

  virtual const HierarchicalMesh*& GetHierarchicalMeshPointer() { return _HM; }
  virtual const HierarchicalMesh*  GetHierarchicalMesh() const  { return _HM; }

  //
  /// for gmres
  //
  virtual void DeleteVector(VectorInterface& p) const;

  double ScalarProduct(const VectorInterface& y, const VectorInterface& x) const;
  void   Equ(VectorInterface& dst, double s, const VectorInterface& src) const;
  void   Add(VectorInterface& dst, double s, const VectorInterface& src) const;
  void   SAdd(double s1, VectorInterface& dst, double s2, const VectorInterface& src) const;
  double Norm(const VectorInterface& dst) const;  
};
}

#endif



