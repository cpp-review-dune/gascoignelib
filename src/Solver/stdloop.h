#ifndef  __StdLoop_h
#define  __StdLoop_h

#include  "extrapolator.h"
#include  "multilevelsolverinterface.h"
#include  "adaptordata.h"
#include  "meshagentinterface.h"
#include  "monitor.h"
#include  "visualization.h"
#include  "visudatacompvector.h"
#include  "visudatanvector.h"
#include  "functionalmanager.h"

#include  "problemdescriptorinterface.h"
#include  "stdiomanager.h"
#include  "stopwatch.h"

//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

class StdLoop
{
private:

  FunctionalManager*         _FMP;
  MeshAgentInterface*        _MA;
  MultiLevelSolverInterface* _ML;

  const ProblemDescriptorInterface* _PD;

  void WriteMeshAndSolution(const std::string& filename, const NewMultiLevelGhostVector& u) const;
  void WriteSolution(const NewMultiLevelGhostVector& u) const;
  void WriteMesh() const;
  
protected:

  FunctionalManager*& GetFunctionalManagerPointer() { return _FMP;}
  MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  MultiLevelSolverInterface*& GetMultiLevelSolverPointer() { return _ML;}

  const FunctionalManager* GetFunctionalManager() const { return _FMP;}
  const MeshAgentInterface* GetMeshAgent() const { return _MA;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const { return _ML;}

  FunctionalManager* GetFunctionalManager() { return _FMP;}
  MeshAgentInterface* GetMeshAgent() { return _MA;}
  MultiLevelSolverInterface* GetMultiLevelSolver() { return _ML;}


  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_functionals, _clock_write, _clock_estimate;

  int    _niter, _nmax, _nmin, _coarse, _iter;
  double _p;

  std::string _paramfile, _refiner, _estimator, _extrapolate, _reload, _initial;

  Extrapolator      Extra;
  Monitor           Mon;
  StdIoManager      IOM;

  const ProblemDescriptorInterface* GetProblemDescriptor() const { return _PD;}

  nvector<double>    _JErr;
  CompVector<double> _GlobalErr;

  std::vector<std::string> GetAllFunctionalNames () const {
    assert(_FMP);
    return _FMP->GetFunctionalNames();
  }
  std::vector<std::string> GetGridFunctionalNames() const {
    assert(_FMP);
    return _FMP->GetGridFunctionalNames();
  }

  virtual void NewFunctionalManager();


  // new vectors

  virtual std::string Solve(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, string name="Results/solve");
  virtual void StdLoop::Output(const NewMultiLevelGhostVector& u, string name="Results/solve") const;

  nvector<double> ComputeAllFunctionals(NewMultiLevelGhostVector& f, NewMultiLevelGhostVector& u) const;
  nvector<double> GetExactValues() const;

  virtual void ComputeGlobalErrors(const NewMultiLevelGhostVector& u);
  virtual void EtaVisu(std::string name, int i, const nvector<double>& eta);
  virtual void AdaptMesh(const nvector<double>& eta);
  virtual nvector<double> Functionals(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f);
  virtual double Estimator(nvector<double>& eta, NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f);

  virtual void InitSolution(NewMultiLevelGhostVector& u);
  virtual void CopyVector(GlobalVector& dst, NewMultiLevelGhostVector& src);
  virtual void CopyVector(NewMultiLevelGhostVector& dst, GlobalVector& src);

public:

  StdLoop();
  virtual ~StdLoop();

  virtual void BasicInit(const std::string& pfile, const ProblemDescriptorInterface& PD);

  void run();
};

/*-----------------------------------------*/

#endif
