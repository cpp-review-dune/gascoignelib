#ifndef  __StdMultiLevelSolver_h
#define  __StdMultiLevelSolver_h
//////////////////////////////////////////////
///
///@brief
/// Default nonlinear MultilevelSolver

/// - stores MultiGridMeshInterace
/// - stores array of MGInterpolator
/// - stores array of SolverInterface
///
//////////////////////////////////////////////

#include  "multilevelsolverinterface.h"
#include  "multilevelsolverdata.h"
#include  "problemdescriptorinterface.h"
#include  "monitor.h"
#include  "stopwatch.h"
#include  "mginterpolatorinterface.h"
#include  "multilevelghostvector.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
class StdMultiLevelSolver : public MultiLevelSolverInterface
{
  private :

  std::vector<SolverInterface*>  _SP;
  const MeshAgentInterface* _MAP;
  std::vector<MgInterpolatorInterface*>   _Interpolator;

  protected :

  const MeshAgentInterface* GetMeshAgent() const {return _MAP;}

  mutable MultiLevelGhostVector _cor, _res, _mg0, _mg1;
  std::set<MultiLevelGhostVector>  _MlVectors;

  mutable StopWatch   _clock_residual, _clock_solve;

  mutable int ComputeLevel;
  mutable int oldnlevels;

  const ParamFile*  _paramfile;

  Monitor*         MON;
  MultiLevelSolverData*          DataP;
  const ProblemDescriptorInterface*      _PD;

  virtual void NewSolvers();

  virtual SolverInterface* NewSolver(int solverlevel);
  virtual void NewMgInterpolator();
  virtual void SolverNewMesh();

  virtual void SetComputeLevel(int level) {ComputeLevel=level;}

  virtual double NewtonNorm(const MultiLevelGhostVector& u) const {
    return GetSolver(ComputeLevel)->NewtonNorm(u(ComputeLevel));
  }
  virtual void mgstep(std::vector<double>& res, std::vector<double>& rw, int l, int maxl, int minl, std::string& p0, std::string p, MultiLevelGhostVector& u, MultiLevelGhostVector& b, MultiLevelGhostVector& v);

  virtual void Cg   (MultiLevelGhostVector& x, const MultiLevelGhostVector& f, CGInfo& info);
  virtual void Gmres(MultiLevelGhostVector& x, const MultiLevelGhostVector& f, CGInfo& info);

  virtual void ViewProtocoll() const;

 public:

  // Constructor

  StdMultiLevelSolver();
  ~StdMultiLevelSolver();

  std::string GetName() const {return "StdMultiLevelSolver";}

  void RegisterVectorAndMemory(const MultiLevelGhostVector& g);
  void RegisterVectorAndMemory();
  void RegisterVector(MultiLevelGhostVector& g);
  void RegisterMatrix();
  void ReInitMatrix();
  void ReInitVector();

  void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile);

  // Zugriff

  virtual void SetState(const std::string& s) {
    for(int l=0;l<_SP.size();l++) _SP[l]->SetState(s);
  }

//  const CGInfo& GetLinearInfo   () const { return DataP->GetLInfo();}
//  const NLInfo& GetNonlinearInfo() const { return DataP->GetNLInfo();}
  const ProblemDescriptorInterface* GetProblemDescriptor() const { return _PD;}

  int nlevels()                 const { assert(GetMeshAgent()); return GetMeshAgent()->nlevels();}
  virtual int FinestLevel  ()  const { return nlevels()-1;}
  virtual int CoarsestLevel()  const { return 0;}

  SolverInterface*& GetSolverPointer(int l) {assert(l<_SP.size()); return _SP[l];}
  SolverInterface* GetSolver(int l) {assert(l<_SP.size()); return _SP[l];}
  const SolverInterface* GetSolver(int l) const {assert(l<_SP.size()); return _SP[l];}
  SolverInterface* GetSolver() {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}
  const SolverInterface* GetSolver() const {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}

  void SetMonitorPtr(Monitor* mon) { MON = mon;}

  void ReInit(const ProblemDescriptorInterface& PDX);
  void SetProblem(const ProblemDescriptorInterface& PDX);

  // neue vektoren

  std::string Solve(int level, MultiLevelGhostVector& x, const MultiLevelGhostVector& b);
  void InterpolateSolution(MultiLevelGhostVector& u, const GlobalVector& uold) const;

  virtual void NewtonVectorZero(MultiLevelGhostVector& w) const;
  virtual double NewtonResidual(MultiLevelGhostVector& y, const MultiLevelGhostVector& x, const MultiLevelGhostVector& b) const;
  virtual double NewtonUpdate(double& rr, MultiLevelGhostVector& x, MultiLevelGhostVector& dx, MultiLevelGhostVector& r, const MultiLevelGhostVector& f, NLInfo& nlinfo);
  virtual void NewtonLinearSolve(MultiLevelGhostVector& x, const MultiLevelGhostVector& b, const MultiLevelGhostVector& u, CGInfo& info);
  virtual void NewtonMatrixControl(MultiLevelGhostVector& u, NLInfo& nlinfo);

  virtual void AssembleMatrix(MultiLevelGhostVector& u, NLInfo& nlinfo);
  void AssembleMatrix(MultiLevelGhostVector& u);
  /// not used in the library -- might be used in local
  virtual void ComputeIlu(MultiLevelGhostVector& u);
  
  virtual void BoundaryInit(MultiLevelGhostVector& u) const;

  virtual void SolutionTransfer(int high, int low, MultiLevelGhostVector& u) const;
  virtual void Transfer(int high, int low, MultiLevelGhostVector& u) const;
  virtual void SolutionTransfer(MultiLevelGhostVector& u) const;
  
  void precondition(MultiLevelGhostVector& x, MultiLevelGhostVector& y);
  void vmulteqgmres(MultiLevelGhostVector& y, const MultiLevelGhostVector&  x) const;
  
  virtual void LinearMg(int minlevel, int maxlevel, MultiLevelGhostVector& u, const MultiLevelGhostVector& f, CGInfo&);

  virtual void NewtonOutput(NLInfo& nlinfo) const;

  double ComputeFunctional(MultiLevelGhostVector& f, const MultiLevelGhostVector& u, const Functional* FP) const;

  void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const;
  void Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const;
};
}

#endif


