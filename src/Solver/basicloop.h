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
#include  "multilevelghostvector.h"
#include  "solverinfos.h"

//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

namespace Gascoigne
{
class BasicLoop
{
private:

  MeshAgentInterface*        _MA;
  MultiLevelSolverInterface* _ML;
  SolverInfos*               _SI;
  
  void WriteMeshAndSolution(const std::string& filename, const MultiLevelGhostVector& u) const;
  void WriteSolution(const MultiLevelGhostVector& u) const;
  void WriteMesh() const;
  
protected:

  MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  MultiLevelSolverInterface*& GetMultiLevelSolverPointer() { return _ML;}

  const MeshAgentInterface* GetMeshAgent() const { return _MA;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const { return _ML;}

  MeshAgentInterface* GetMeshAgent() { return _MA;}
  MultiLevelSolverInterface* GetMultiLevelSolver() { return _ML;}

        SolverInfos*& GetSolverInfosPointer() { return _SI;}
        SolverInfos* GetSolverInfos()         { return _SI;}
  const SolverInfos* GetSolverInfos() const   { return _SI;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_functionals, _clock_write, _clock_estimate;

  int    _niter, _nmax, _nmin, _coarse, _iter;
  double _p;

  std::string _reload, _initial;
  const ParamFile*  _paramfile;
  

  Monitor           Mon;
  StdIoManager      IOM;

  GlobalVector _GlobalErr;

  // new vectors

  virtual std::string Solve(MultiLevelGhostVector& u, MultiLevelGhostVector& f, std::string name="Results/solve");

  virtual void PrintMeshInformation(int outputlevel=0) const;
  virtual void Output(const MultiLevelGhostVector& u, std::string name="Results/solve") const;

  virtual void ComputeGlobalErrors(const MultiLevelGhostVector& u);

  virtual void InitSolution(MultiLevelGhostVector& u);
  virtual void CopyVector(GlobalVector& dst, MultiLevelGhostVector& src);
  virtual void CopyVector(MultiLevelGhostVector& dst, GlobalVector& src);

public:

  BasicLoop();
  virtual ~BasicLoop();

  virtual void BasicInit(const ParamFile* paramfile);

  void run(const ProblemDescriptorInterface* PD);
  void ClockOutput() const;
};
}

/*-----------------------------------------*/

#endif
