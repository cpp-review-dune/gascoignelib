#ifndef  __MultiLevelSolverInterface_h
#define  __MultiLevelSolverInterface_h


#include  "meshagentinterface.h"
#include  "solverinterface.h"
#include  "monitor.h"
#include  "paramfile.h"
#include  "nlinfo.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments MultiLevelSolverInterface

  ///
  ///
  /////////////////////////////////////////////

  class MultiLevelGhostVector;

  class MultiLevelSolverInterface
  {
    private:

    protected:

    public:
      MultiLevelSolverInterface() {}
      virtual ~MultiLevelSolverInterface() {}

      virtual std::string GetName() const=0;
      virtual void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile)=0;
      virtual void SetProblem(const ProblemDescriptorInterface& PDX)=0;
      virtual void ReInit(const ProblemDescriptorInterface& PDX)=0;
      virtual void SetMonitorPtr(Monitor* mon)=0;

      virtual void ReInitMatrix()=0;
      virtual void ReInitVector()=0;

      virtual int nlevels() const=0;

      virtual SolverInterface* GetSolver(int l)=0;
      virtual const SolverInterface* GetSolver(int l) const=0;
      virtual SolverInterface* GetSolver()=0;
      virtual const SolverInterface* GetSolver() const=0;

//      virtual void SetState(const std::string& s)=0;
      virtual void AssembleMatrix(MultiLevelGhostVector& u, NLInfo& nlinfo)=0;
      virtual void AssembleMatrix(MultiLevelGhostVector& u)=0;
      virtual void ComputeIlu(MultiLevelGhostVector& u)=0;
      virtual void ComputeIlu()=0;
      
      virtual void BoundaryInit(MultiLevelGhostVector& u) const=0;

      //
      /// vector - manamgement
      //

      virtual void RegisterVector(MultiLevelGhostVector& g)=0;
      virtual void MemoryVector(MultiLevelGhostVector& g)=0;
      virtual void DeleteVector(MultiLevelGhostVector& g)=0;
      //virtual void RegisterVectorAndMemory(const MultiLevelGhostVector& g)=0;
      virtual void RegisterVectorsOnSolvers()=0;
      virtual void RegisterMatrix()=0;

      //
      /// vector 
      //

      virtual std::string LinearSolve(int level, MultiLevelGhostVector& u, const MultiLevelGhostVector& b, CGInfo& info)=0;
      virtual std::string LinearSolve(MultiLevelGhostVector& u, const MultiLevelGhostVector& b, CGInfo& info) {
        return LinearSolve(nlevels()-1,u,b,info);
      }
      virtual std::string Solve(int level, MultiLevelGhostVector& x, const MultiLevelGhostVector& b, NLInfo& nlinfo)=0;
      virtual std::string Solve(MultiLevelGhostVector& x, const MultiLevelGhostVector& b, NLInfo& nlinfo) {
        return Solve(nlevels()-1,x,b,nlinfo);
      }
      virtual void InterpolateSolution(MultiLevelGhostVector& u, const GlobalVector& uold) const=0;
      virtual double ComputeFunctional(MultiLevelGhostVector& f, const MultiLevelGhostVector& u, const Functional* FP) const=0;
      virtual void Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const=0;
      virtual void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const=0;
      virtual void AssembleDualMatrix(MultiLevelGhostVector& u)=0;
      virtual void vmulteq(MultiLevelGhostVector& y, const MultiLevelGhostVector&  x) const=0;
  };
}

#endif
