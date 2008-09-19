#ifndef __cg_h 
#define __cg_h 

#include "cginfo.h" 

/********************************************************************/

// S = Solver
// V = Vector

namespace Gascoigne
{
template <class S, class V>
class CG
{
  S   &solver;

public:

  CG(S& s) : solver(s) {}

  void solve(V& x, const V& b, CGInfo& info);
};
}

/********************************************************************/


#endif 
