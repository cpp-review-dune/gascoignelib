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

#include  "problemdescriptorinterface.h"
#include  "stdiomanager.h"
#include  "stopwatch.h"
#include  "paramfile.h"

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

  MeshAgentInterface*        _MA;
  MultiLevelSolverInterface* _ML;

  vector<const Functional*>   _FV;

  void WriteMeshAndSolution(const std::string& filename, const NewMultiLevelGhostVector& u) const;
  void WriteSolution(const NewMultiLevelGhostVector& u) const;
  void WriteMesh() const;
  
protected:

  MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  MultiLevelSolverInterface*& GetMultiLevelSolverPointer() { return _ML;}

  const vector<const Functional*>& GetFunctionals() const { return _FV;}
  const MeshAgentInterface* GetMeshAgent() const { return _MA;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const { return _ML;}

  MeshAgentInterface* GetMeshAgent() { return _MA;}
  MultiLevelSolverInterface* GetMultiLevelSolver() { return _ML;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_functionals, _clock_write, _clock_estimate;

  int    _niter, _nmax, _nmin, _coarse, _iter;
  double _p;

  const ParamFile*  _paramfile;
  std::string _refiner, _estimator, _extrapolate, _reload, _initial;

  Extrapolator      Extra;
  Monitor           Mon;
  StdIoManager      IOM;

  nvector<double>    _JErr;
  CompVector<double> _GlobalErr;

  // new vectors

  virtual std::string Solve(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, string name="Results/solve");
  virtual void StdLoop::Output(const NewMultiLevelGhostVector& u, string name="Results/solve") const;

  nvector<double> ComputeFunctionals(NewMultiLevelGhostVector& f, NewMultiLevelGhostVector& u, const vector<const Functional*>& J) const;

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

  virtual void BasicInit(const ParamFile* paramfile);

  void SetFunctionals(const vector<const Functional*>& fv) { _FV =  fv;}

  void run(const ProblemDescriptorInterface* PD);
};

/*-----------------------------------------*/

#endif
