#ifndef  __NumericInterface_h
#define  __NumericInterface_h

#include <iostream>
#include  <cstdlib>

namespace Gascoigne
{

class DiscretizationInterface;
class SolverInterface;
class MeshAgentInterface;

/*----------------------------------------------------------------------------*/

class NumericInterface
{
public:
  NumericInterface() { };
  virtual ~NumericInterface() { };

  virtual DiscretizationInterface* NewDiscretization(int level=0) const
  { std::cerr << "NumericInterface incomplete: NewDiscretization" << std::endl; abort();}

  virtual SolverInterface* NewSolver(int level=0) const
  { std::cerr << "NumericInterface incomplete: NewSolver" << std::endl; abort(); }

  virtual MeshAgentInterface* NewMeshAgent() const
  { std::cerr << "NumericInterface incomplete: NewMeshAgent" << std::endl; abort(); }
};

}

#endif
