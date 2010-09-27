#ifndef  __MultiLevelAlgorithm_h
#define  __MultiLevelAlgorithm_h

#include  "algorithm.h"
#include  "multilevelsolver.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class MultiLevelAlgorithm : public Algorithm
{
 private:

  MultiLevelSolver*           _S;

  int                         _coarselevel;

  std::string                 _mgtype;
  double                      _mgomega;
  
 protected:

  virtual const SolverInterface* GetSolver() const { return _S->GetSolver();}
  virtual       SolverInterface* GetSolver()       { return _S->GetSolver();}

  virtual const SolverInterface* GetSolver(int i) const { return _S->GetSolver(i);}
  virtual       SolverInterface* GetSolver(int i)       { return _S->GetSolver(i);}

  virtual const MultiLevelSolver* GetMultiLevelSolver() const { return _S;}
  virtual       MultiLevelSolver* GetMultiLevelSolver()       { return _S;}

  void  DeleteVector(VectorInterface& u) const { _S->DeleteVector(u);}
  void  AssembleMatrixAndIlu(VectorInterface& u);
  void  LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo);
  void  NonLinear(VectorInterface& u, VectorInterface& f,
		  const std::string& problemlabel, int iter);
  void VWCycle(std::vector<double>& res, std::vector<double>& rw, 
	       int l, int finelevel, int coarselevel, const std::string& p,
	       VectorInterface& u, VectorInterface& b, VectorInterface& v);
  void LinearMGSolve     (VectorInterface& du, const VectorInterface& y, CGInfo& cginfo);
  virtual void Precondition(VectorInterface& x, VectorInterface& y);

public:

  MultiLevelAlgorithm() :  _S(NULL) {}
  virtual ~MultiLevelAlgorithm() {}

  void BasicInit(const ParamFile* paramfile, MultiLevelSolver* MLS, const NumericInterface* NI,
		 const ProblemContainer* PC);

  void ReInitVector(VectorInterface& u) const { _S->ReInitVector(u);} 
  void RunLinear(const std::string& problemlabel);
  void RunNonLinear(const std::string& problemlabel,int iter=0);
  void GlobalRefineLoop(const std::string& problemlabel);
  void LocalRefineLoop(const std::string& problemlabel, FunctionalContainer* FC=NULL);
};
}

/*-----------------------------------------*/

#endif
