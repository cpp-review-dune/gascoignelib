#ifndef  __MultiLevelSolverInterface_h
#define  __MultiLevelSolverInterface_h


#include  "meshagentinterface.h"
#include  "solverinterface.h"
#include  "monitor.h"
#include  "paramfile.h"
#include  "nlinfo.h"
#include  "vectorinterface.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments MultiLevelSolverInterface

  ///
  ///
  /////////////////////////////////////////////

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
      virtual void AssembleMatrix(VectorInterface& u, NLInfo& nlinfo)=0;
      virtual void AssembleMatrix(VectorInterface& u)=0;
      virtual void ComputeIlu(VectorInterface& u)=0;
      virtual void ComputeIlu()=0;
      
      virtual void BoundaryInit(VectorInterface& u) const=0;

      //
      /// vector - manamgement
      //

      virtual void RegisterVector(VectorInterface& g)=0;
      virtual void RegisterVectorsOnSolvers()=0;
      virtual void RegisterMatrix()=0;

      //
      /// vector 
      //

      virtual std::string LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info)=0;
      virtual std::string LinearSolve(VectorInterface& u, const VectorInterface& b, CGInfo& info) {
        return LinearSolve(nlevels()-1,u,b,info);
      }
      virtual std::string Solve(int level, VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo)=0;
      virtual std::string Solve(VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo) {
        return Solve(nlevels()-1,x,b,nlinfo);
      }
      virtual void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const=0;
      virtual double ComputeFunctional(VectorInterface& f, const VectorInterface& u, const Functional* FP) const=0;
      virtual void Transfer(int l, GlobalVector& ul, const GlobalVector& uf) const=0;
      virtual void SolutionTransfer(int l, GlobalVector& ul, const GlobalVector& uf) const=0;
      virtual void AssembleDualMatrix(VectorInterface& u)=0;
      virtual void vmulteq(VectorInterface& y, const VectorInterface&  x) const=0;

      virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src)const=0;
      virtual void Zero(VectorInterface& dst)const=0;
  };
}

#endif
