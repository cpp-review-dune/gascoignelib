#ifndef  __StdMultiLevelSolver_h
#define  __StdMultiLevelSolver_h

#include  "multilevelsolverinterface.h"
#include  "stdmultilevelsolverdata.h"
#include  "problemdescriptorinterface.h"
#include  "monitor.h"
#include  "stopwatch.h"
#include  "mginterpolatorinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear MultilevelSolver

/// - stores MultiGridMeshInterace
/// - stores array of MGInterpolator
/// - stores array of SolverInterface
///
//////////////////////////////////////////////

class StdMultiLevelSolver : public MultiLevelSolverInterface
{
  private :

  std::vector<SolverInterface*>  _SP;
  const MeshAgentInterface* _MAP;
  std::vector<MgInterpolatorInterface*>   _Interpolator;

  protected :

  const MeshAgentInterface* GetMeshAgent() const {return _MAP;}
  std::vector<SolverInterface*>& GetSolverPointers() { return _SP; }
  const std::vector<SolverInterface*>& GetSolverPointers() const { return _SP; }


  mutable VectorInterface _cor, _res, _mg0, _mg1;

  mutable StopWatch   _clock_residual, _clock_solve;

  mutable int ComputeLevel;
  mutable int oldnlevels;

  const ParamFile*  _paramfile;

  Monitor*                          MON;
  StdMultiLevelSolverData*          DataP;
  const ProblemDescriptorInterface*      _PD;

  virtual void NewSolvers();

  virtual SolverInterface* NewSolver(int solverlevel);
  virtual void NewMgInterpolator();
  virtual void SolverNewMesh();

  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const { return _PD;}
  virtual SolverInterface*& GetSolverPointer(int l) {assert(l<_SP.size()); return _SP[l];}
  virtual void SetComputeLevel(int level) {ComputeLevel=level;}

  virtual double NewtonNorm(const VectorInterface& u) const {
    return GetSolver(ComputeLevel)->NewtonNorm(u);
  }
  virtual void mgstep(std::vector<double>& res, std::vector<double>& rw, int l, int maxl, int minl, std::string& p0, std::string p, VectorInterface& u, VectorInterface& b, VectorInterface& v);

  virtual void Cg   (VectorInterface& x, const VectorInterface& f, CGInfo& info);
  virtual void Gmres(VectorInterface& x, const VectorInterface& f, CGInfo& info);

  virtual void ViewProtocoll() const;

  virtual void SolutionTransfer(VectorInterface& u) const;
  virtual void SolutionTransfer(int high, int low, VectorInterface& u) const;
  virtual void Transfer(int high, int low, VectorInterface& u) const;

 public:

  // Constructor

  StdMultiLevelSolver();
  ~StdMultiLevelSolver();

  std::string GetName() const {return "StdMultiLevelSolver";}

  void RegisterVectors();
  void RegisterMatrix();
  void ReInitMatrix();
  void ReInitVectors();
  void ReInitVector(VectorInterface& v);

  void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile);

  // Zugriff

//  virtual void SetState(const std::string& s) {
//    for(int l=0;l<_SP.size();l++) _SP[l]->SetState(s);
//  }


  int nlevels()                 const { assert(GetMeshAgent()); return GetMeshAgent()->nlevels();}
  virtual int FinestLevel  ()  const { return nlevels()-1;}
  virtual int CoarsestLevel()  const { return 0;}

  
  
  SolverInterface* GetSolver(int l) {assert(l<_SP.size()); return _SP[l];}
  const SolverInterface* GetSolver(int l) const {assert(l<_SP.size()); return _SP[l];}
  SolverInterface* GetSolver() {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}
  const SolverInterface* GetSolver() const {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}

  void SetMonitorPtr(Monitor* mon) { MON = mon;}

  void ReInit(const ProblemDescriptorInterface& PDX);
  void SetProblem(const ProblemDescriptorInterface& PDX);

  // neue vektoren

  std::string LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info);
  std::string Solve(int level, VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo);
  void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const;

  virtual void NewtonVectorZero(VectorInterface& w) const;
  virtual double NewtonResidual(VectorInterface& y, const VectorInterface& x, const VectorInterface& b) const;
  virtual double NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo);
/*   virtual void NewtonUpdateShowCompResiduals(VectorInterface& x, VectorInterface& r, const VectorInterface& f); */
  virtual void NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info);
  virtual void NewtonMatrixControl(VectorInterface& u, NLInfo& nlinfo);
  virtual void NewtonOutput(NLInfo& nlinfo) const;

  void AssembleMatrix(VectorInterface& u, NLInfo& nlinfo);
  void AssembleMatrix(VectorInterface& u);
  /// not used in the library -- might be used in local
  void ComputeIlu(VectorInterface& u);
  void ComputeIlu();
  
  void BoundaryInit(VectorInterface& u) const;
  
  void vmulteq(VectorInterface& y, const VectorInterface&  x) const;
  
  virtual void LinearMg(int minlevel, int maxlevel, VectorInterface& u, const VectorInterface& f, CGInfo&);

  double ComputeFunctional(VectorInterface& f, const VectorInterface& u, const Functional* FP) const;

  void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const;
  void Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const;
  void AssembleDualMatrix(VectorInterface& u);

  // fuer gmres
  
  virtual void precondition(VectorInterface& x, VectorInterface& y);
  virtual void DeleteVector(VectorInterface& p);
  virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src)const;
  void Zero(VectorInterface& dst)const;

  void AddNodeVector(const std::string& name, VectorInterface& q);
  void DeleteNodeVector(const std::string& q);

  void newton(VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info);
};
}

#endif


