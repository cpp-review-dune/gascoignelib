#ifndef  __LpsEquation_h
#define  __LpsEquation_h

#include "equation.h"
#include "stabinterface.h"


namespace Gascoigne
{

/*-----------------------------------------*/

  //////////////////////////////////////////////
  //
  ///@brief
  /// Interface class for Lps Elements
  ///
  ///
  //////////////////////////////////////////////

  class LpsEquation : public virtual Equation
  {
    public:
      LpsEquation() {}
      ~LpsEquation() {}

      virtual StabInterface& GetStabilization() const { assert(0);}

      virtual void init(const nmatrix<double>& H, const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
      virtual void init(const nmatrix<double>& H, const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
   
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

  virtual void StabilizationResidual(LocalVector& F, const FemFunction& U, const FemFunction& UP, const FemFunction& N, const FemFunction& NP) const
    {
      assert(0);
    }

  virtual void StabilizationMatrix(EntryMatrix& A, const FemFunction& U, const FemFunction& UP, const FemFunction& M, const FemFunction& MP, const FemFunction& N, const FemFunction& NP) const
    {
      assert(0);
    };
  };
}

#endif
