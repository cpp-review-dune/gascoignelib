#ifndef  __LpsEquation_h
#define  __LpsEquation_h

#include "equation.h"


namespace Gascoigne
{
  
  //////////////////////////////////////////////
  //
  ///@brief
  /// Interface class for Lps Elements
  ///
  ///
  //////////////////////////////////////////////

  class LpsEquation : public virtual Equation
  {
    private:

    protected:

    public:
      LpsEquation() {}
      ~LpsEquation() {}

      virtual void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
      virtual void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
   
      virtual void lpspointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
        lpspoint(h,U,v);
      }
      virtual void lpspointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
        lpspoint(h,U,v);
      }

      virtual void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const=0;
      virtual void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const DerivativeVector& Mp) const=0;
  };
}

#endif
