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

  void NewMeshInterpretor(int dimension, const std::string& discname);

public:


//
////  Con(De)structor 
//

  LocalSolver() : StdSolver() {}
  ~LocalSolver() {}

  std::string GetName() const {return "Local";}

  int n() const {return GetMeshInterpretor()->n();}

  void InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const;
  void HNZero(GlobalVector& x) const;
  void HNAverage(GlobalVector& x) const;

  void AddGlobalData(const GlobalVector& d);
};


#endif
