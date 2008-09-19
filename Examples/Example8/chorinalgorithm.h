#ifndef  __ChorinAlgorithm_h
#define  __ChorinAlgorithm_h

#include  "splittingsolver.h"
#include  "multilevelsolver.h"
#include  "nonstationaryalgorithm.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class ChorinAlgorithm : public NonstationaryAlgorithm
{
 protected:

  int         niter;
  std::string initial;

  const SplittingSolver* GetSplittingSolver() const 
    { 
      const SplittingSolver* S = dynamic_cast<const SplittingSolver*>(GetSolver());
      assert(S);
      return S;
    }
       SplittingSolver* GetSplittingSolver()
    { 
      SplittingSolver* S = dynamic_cast<SplittingSolver*>(GetSolver());
      assert(S);
      return S;
    }
  void ChorinUzawa(const std::string& problem1, const std::string& problem2);
  void VanKan(const std::string& problem1, const std::string& problem2);

public:

           ChorinAlgorithm() : NonstationaryAlgorithm(), niter(0),initial("")  {}
  virtual ~ChorinAlgorithm() {}

  void Run(const std::string& problem1, const std::string& problem2);
};
}

/*-----------------------------------------*/

#endif
