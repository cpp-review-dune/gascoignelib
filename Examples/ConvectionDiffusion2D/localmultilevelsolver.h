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

 SolverInterface* NewSolver(int solverlevel);

public:


//
////  Con(De)structor 
//

  LocalMultiLevelSolver() : StdMultiLevelSolver() {}
 ~LocalMultiLevelSolver() {}

  std::string GetName() const {return "Local";}

  void SolutionTransfer(int level, GlobalVector& uc, GlobalVector& uf) const;

};


#endif
