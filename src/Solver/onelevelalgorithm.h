#ifndef  __OneLevelAlgorithm_h
#define  __OneLevelAlgorithm_h

#include  "algorithm.h"
#include  "stdsolver.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class OneLevelAlgorithm : public Algorithm
{
 private:

  SolverInterface*        _S;
  const ProblemContainer* _PC;
  
 protected:

  virtual const SolverInterface* GetSolver() const { return _S;}
  virtual       SolverInterface* GetSolver()       { return _S;}

  void IluSolver(VectorInterface& du, const VectorInterface& f, CGInfo& info);

  void  ReInitVector(VectorInterface& u) const { _S->ReInitVector(u);} 
  void  DeleteVector(VectorInterface& u) const { _S->DeleteVector(u);}
  void  AssembleMatrixAndIlu(VectorInterface& u);
  void  LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo);

public:

  OneLevelAlgorithm() :  _S(NULL) {}
  virtual ~OneLevelAlgorithm() { if (_S) delete _S;}

  virtual void BasicInit(const ParamFile* paramfile, const NumericInterface* NI,
			 const ProblemContainer* PC);

  void RunLinear   (const std::string& problemlabel);
  void RunNonLinear(const std::string& problemlabel);
};
}

/*-----------------------------------------*/

#endif
