#ifndef __NewmarkSolver_h
#define __NewmarkSolver_h

#include  "timesolver.h"

/*----------------------------------------------------------------------------*/

class NewmarkSolver : public Gascoigne::TimeSolver
{
public:

  NewmarkSolver() : Gascoigne::TimeSolver() {}
  void Form(Gascoigne::VectorInterface& gy, const Gascoigne::VectorInterface& gx, double d) const;
  void FormWithoutMass(Gascoigne::VectorInterface& gy, const Gascoigne::VectorInterface& gx, 
		       double d, double s=1.) const;
  void AssembleMatrix(const Gascoigne::VectorInterface& gu, double d);
};

/*----------------------------------------------------------------------------*/

#endif
