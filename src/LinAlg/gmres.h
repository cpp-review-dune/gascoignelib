#ifndef __gmresclass_h
#define __gmresclass_h

#include "mult.h"
#include "cginfo.h"
#include "nvector.h"

/********************************************************************/

namespace Gascoigne
{
template<class SOLVER, class PRECONDITIONER, class VECTOR>
class GMRES
{
  typedef nvector<double>  dvector;

  nvector<dvector> H;
  nvector<double>  gamma, ci, si;

  std::vector<VECTOR> mem;
  
  int vmax, left_precondition;
  
  void   new_memory          ();
  void   givens_rotation     (dvector&, int);
  void   solution            (VECTOR&, VECTOR&, int);
  double orthogonalization   (dvector&, int, VECTOR&) const;
  bool   reortho_test        (const VECTOR&, double) const;

  SOLVER&         A;
  PRECONDITIONER& P;
  
 public:

  GMRES(SOLVER&, PRECONDITIONER&, int);
  ~GMRES();
  void init();

  int solve          (VECTOR& x, const VECTOR& b, CGInfo& info);
  int restarted_solve(VECTOR& x, const VECTOR& b, CGInfo& info);
};
}

/********************************************************************/

#endif
