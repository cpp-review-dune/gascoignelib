#ifndef  __LocalSolver_h
#define  __LocalSolver_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalSolver

////
////
/////////////////////////////////////////////


#include  "stdsolver.h"

class LocalSolver : public StdSolver
{
private:


protected:


public:


//
////  Con(De)structor 
//

  LocalSolver();
  ~LocalSolver() {}

  MeshInterpretorInterface* NewMeshInterpretor(int dimension, const std::string& discname);

};


#endif
