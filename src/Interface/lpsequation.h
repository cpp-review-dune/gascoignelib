#ifndef  __LpsEquation_h
#define  __LpsEquation_h

#include "equation.h"

//////////////////////////////////////////////
//
///@brief
/// Interface class for Lps Elements
///
///
//////////////////////////////////////////////

namespace Gascoigne
{

/*-----------------------------------------*/

class LpsEquation : public virtual Equation
{
 protected:

public:

  LpsEquation() {}
  ~LpsEquation() {}

  virtual void lpspoint
    (double h, const FemFunction& U, FemData& Q, const Vertex2d& v) const { assert(0);}
  virtual void lpspoint
    (double h, const FemFunction& U, FemData& Q, const Vertex3d& v) const { assert(0);}
 
  virtual void lpspointmatrix
    (double h, const FemFunction& U, FemData& Q, const Vertex2d& v) const
    {
      lpspoint(h,U,Q,v);
    };
  virtual void lpspointmatrix
    (double h, const FemFunction& U, FemData& Q, const Vertex3d& v) const
    {
      lpspoint(h,U,Q,v);
    };

  virtual void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const=0;

  virtual void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const DerivativeVector& Mp) const=0;
};

}

#endif
