#ifndef  __BasicLoop_h
#define  __BasicLoop_h

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

class BasicLoop
{
private:

  MeshAgentInterface*        _MA;
  MultiLevelSolverInterface* _ML;

  void WriteMeshAndSolution(const std::string& filename, const NewMultiLevelGhostVector& u) const;
  void WriteSolution(const NewMultiLevelGhostVector& u) const;
  void WriteMesh() const;
  
protected:

  MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  MultiLevelSolverInterface*& GetMultiLevelSolverPointer() { return _ML;}

  const MeshAgentInterface* GetMeshAgent() const { return _MA;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const { return _ML;}

  MeshAgentInterface* GetMeshAgent() { return _MA;}
  MultiLevelSolverInterface* GetMultiLevelSolver() { return _ML;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_functionals, _clock_write, _clock_estimate;

  int    _niter, _nmax, _nmin, _coarse, _iter;
  double _p;

  std::string _reload, _initial;
  const ParamFile*  _paramfile;
  

  Monitor           Mon;
  StdIoManager      IOM;

  CompVector<double> _GlobalErr;

  // new vectors

  virtual std::string Solve(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, std::string name="Results/solve");
  virtual void Output(const NewMultiLevelGhostVector& u, std::string name="Results/solve") const;

  virtual void ComputeGlobalErrors(const NewMultiLevelGhostVector& u);

  virtual void InitSolution(NewMultiLevelGhostVector& u);
  virtual void CopyVector(GlobalVector& dst, NewMultiLevelGhostVector& src);
  virtual void CopyVector(NewMultiLevelGhostVector& dst, GlobalVector& src);

public:

  BasicLoop();
  virtual ~BasicLoop();

  virtual void BasicInit(const Gascoigne::ParamFile* paramfile);

  void run(const ProblemDescriptorInterface* PD);
  void ClockOutput() const;
};

/*-----------------------------------------*/

#endif
