#ifndef  __NumericInterface_h
#define  __NumericInterface_h

#include  "meshagentinterface.h"
#include  "discretizationinterface.h"
#include  "solverinterface.h"

namespace Gascoigne
{

/*----------------------------------------------------------------------------*/

class NumericInterface
{
public:

  virtual ~NumericInterface() {};

  virtual DiscretizationInterface* NewDiscretization(int level=0) const 
  { std::cerr << "NumericInterface incomplete (1)" << std::endl; abort();}

  virtual SolverInterface*         NewSolver(int level=0)     const 
  { std::cerr << "NumericInterface incomplete (2)" << std::endl; abort(); }

  virtual MeshAgentInterface* NewMeshAgent() const
  { std::cerr << "NumericInterface incomplete (4)" << std::endl; abort(); }
};

}

#endif
