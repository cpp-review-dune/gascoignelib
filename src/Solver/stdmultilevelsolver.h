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
#include  "newmultilevelghostvector.h"


using namespace std;

/*-------------------------------------------------------------*/

class StdMultiLevelSolver : public MultiLevelSolverInterface
{
  private :

  vector<SolverInterface*>  _SP;
  const MeshAgentInterface* _MAP;
  vector<MgInterpolatorInterface*>   _Interpolator;

  protected :

  const MeshAgentInterface* GetMeshAgent() const {return _MAP;}

  mutable NewMultiLevelGhostVector _cor, _res, _mg0, _mg1;
  set<NewMultiLevelGhostVector>  _MlVectors;

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

  virtual double NewtonNorm(const NewMultiLevelGhostVector& u) const {
    return GetSolver(ComputeLevel)->NewtonNorm(u(ComputeLevel));
  }
  virtual void mgstep(vector<double>& res, vector<double>& rw, int l, int maxl, int minl, string& p0, string p, NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& b, NewMultiLevelGhostVector& v);

  virtual void Cg   (NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& f, CGInfo& info);
  virtual void Gmres(NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& f, CGInfo& info);

 public:

  // Constructor

  StdMultiLevelSolver();
  ~StdMultiLevelSolver();

  string GetName() const {return "StdMultiLevelSolver";}

  void RegisterVector(NewMultiLevelGhostVector& g) {
//     cerr << "*************registriere:\t"<<g<<endl;
//     _MlVectors.insert(&g);
    _MlVectors.insert(g);
  }
  void MemoryVector();

  void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile);

  // Zugriff

  virtual void SetState(const string& s) {
    for(int l=0;l<_SP.size();l++) _SP[l]->SetState(s);
  }

  const CGInfo& GetLinearInfo   () const { return DataP->info;}
  const CGInfo& GetDualInfo     () const { return DataP->dualinfo;}
  const NLInfo& GetNonlinearInfo() const { return DataP->nlinfo;}
  const ProblemDescriptorInterface* GetProblemDescriptor() const { return _PD;}

  int nlevels()                 const { return GetMeshAgent()->nlevels();}
  virtual int FinestLevel  ()  const { return nlevels()-1;}
  virtual int CoarsestLevel()  const { return 0;}

  SolverInterface* GetSolver(int l) {assert(l<_SP.size()); return _SP[l];}
  const SolverInterface* GetSolver(int l) const {assert(l<_SP.size()); return _SP[l];}
  SolverInterface* GetSolver() {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}
  const SolverInterface* GetSolver() const {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}

  void SetMonitorPtr(Monitor* mon) { MON = mon;}

  void ReInit(const ProblemDescriptorInterface& PDX);
  void NewMesh();
  void SetProblem(const ProblemDescriptorInterface& PDX);

  // neue vektoren

  string Solve(int level, NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& b);
  void InterpolateSolution(NewMultiLevelGhostVector& u, const GlobalVector& uold) const;

  virtual double NewtonResidual(NewMultiLevelGhostVector& y, const NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& b) const;
  virtual double NewtonUpdate(double& rr, NewMultiLevelGhostVector& x, NewMultiLevelGhostVector& dx, NewMultiLevelGhostVector& r, const NewMultiLevelGhostVector& f, NLInfo& nlinfo);
  virtual void NewtonLinearSolve(NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& b, const NewMultiLevelGhostVector& u, CGInfo& info);
  virtual void NewtonMatrixControl(NewMultiLevelGhostVector& u, const NLInfo& nlinfo);

  virtual void AssembleMatrix(NewMultiLevelGhostVector& u);
  /// not used in the library -- might be used in local
  virtual void ComputeIlu(NewMultiLevelGhostVector& u);
  
  virtual void BoundaryInit(NewMultiLevelGhostVector& u) const;

  virtual void SolutionTransfer(int high, int low, NewMultiLevelGhostVector& u) const;
  virtual void SolutionTransfer(NewMultiLevelGhostVector& u) const;
  
  void precondition(NewMultiLevelGhostVector& x, NewMultiLevelGhostVector& y);
  void vmulteqgmres(NewMultiLevelGhostVector& y, const NewMultiLevelGhostVector&  x) const;
  
  virtual void LinearMg(int minlevel, int maxlevel, NewMultiLevelGhostVector& u, const NewMultiLevelGhostVector& f, CGInfo&);

  virtual void NewtonOutput(NLInfo& nlinfo) const;

  double ComputeFunctional(NewMultiLevelGhostVector& f, const NewMultiLevelGhostVector& u, const Functional* FP) const;

  void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const;
};

#endif


