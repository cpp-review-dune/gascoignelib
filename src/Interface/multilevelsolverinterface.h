#ifndef  __MultiLevelSolverInterface_h
#define  __MultiLevelSolverInterface_h


/////////////////////////////////////////////
///
///@brief
///  ... comments MultiLevelSolverInterface

///
///
/////////////////////////////////////////////

#include  "meshagentinterface.h"
#include  "solverinterface.h"
#include  "monitor.h"
//#include  "newmultilevelghostvector.h"

class NewMultiLevelGhostVector;

class MultiLevelSolverInterface
{
public:


private:


protected:


public:

  MultiLevelSolverInterface() {}
  virtual ~MultiLevelSolverInterface() {}

  virtual std::string GetName() const=0;
  virtual void BasicInit(const MeshAgentInterface* GMGM, const std::string& paramfile)=0;
  // temporary
  virtual void BasicInit(MultiGridMeshInterface& GMGM, const std::string& paramfile) {cerr<< "Achtung: Aenderung: jetzt BasicInit(\"const\" MeshAgentInterface& GMGM, const std::string& paramfile)\n"; assert(0);};
  virtual void NewMesh(const ProblemDescriptorInterface* PDX)=0;
  virtual void SetProblem(const ProblemDescriptorInterface& PDX)=0;
  virtual void SetMonitorPtr(Monitor* mon)=0;

  virtual int nlevels() const=0;

  virtual SolverInterface* GetSolver(int l)=0;
  virtual const SolverInterface* GetSolver(int l) const=0;
  virtual SolverInterface* GetSolver()=0;
  virtual const SolverInterface* GetSolver() const=0;

  virtual void SetState(const std::string& s)=0;

  //
  /// vector - manamgement
  //

  virtual void RegisterVector(NewMultiLevelGhostVector& g)=0;

  //
  /// vector 
  //

  virtual std::string Solve(int level, NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& b)=0;
  virtual std::string Solve(NewMultiLevelGhostVector& x, const NewMultiLevelGhostVector& b) {
    return Solve(nlevels()-1,x,b);
  }
  virtual void InterpolateSolution(NewMultiLevelGhostVector& u, const GlobalVector& uold) const=0;
  virtual double ComputeFunctional(NewMultiLevelGhostVector& f, const NewMultiLevelGhostVector& u, const Functional* FP) const=0;
  virtual void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const {assert(0);}
};


#endif
