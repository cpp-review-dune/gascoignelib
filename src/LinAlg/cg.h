#ifndef __cg_h 
#define __cg_h 

#include "cginfo.h" 

/********************************************************************/

// S = Solver
// P = Preconditioner
// V = Vector

template <class S, class P, class V>
class CG
{
  S   &solver;
  P   &precon;

public:

  CG(S& s, P& p) : solver(s), precon(p) {}

  void solve(V& x, const V& b, CGInfo& info);
};

/********************************************************************/


#endif 
