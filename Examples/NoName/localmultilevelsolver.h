#ifndef  __LocalMultiLevelSolver_h
#define  __LocalMultiLevelSolver_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalMultiLevelSolver

////
////
/////////////////////////////////////////////

#include  "stdmultilevelsolver.h"

class LocalMultiLevelSolver : public StdMultiLevelSolver
{
private:


protected:


public:


//
////  Con(De)structor 
//
  LocalMultiLevelSolver();
  ~LocalMultiLevelSolver() {}

  SolverInterface* NewSolver(int solverlevel);

};


#endif
