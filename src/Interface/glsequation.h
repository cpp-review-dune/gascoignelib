#ifndef  __GlsEquation_h
#define  __GlsEquation_h

#include "equation.h"


namespace Gascoigne
{
  
  //////////////////////////////////////////////
  //
  ///@brief
  /// Interface class for  Elements

  ///
  ///
  //////////////////////////////////////////////
  
  class GlsEquation : public virtual Equation
  {
    public:
      GlsEquation() {};
      ~GlsEquation() {};

      //
      /// computation of stabilization parameters
      //
      virtual void glspoint(double h, const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"GlsEquation::glspoint\" not written!" << std::endl;
        abort();
      } 
      
      //
      /// computation of stabilization parameters
      //
      virtual void glspoint(double h, const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"GlsEquation::glspoint\" not written!" << std::endl;
        abort();
      }
     
      //
      /// computation of stabilization parameters for the matrix
      //
      virtual void glspointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
        glspoint(h,U,v);
      }
      
      //
      /// computation of stabilization parameters for the matrix
      //
      virtual void glspointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
        glspoint(h,U,v);
      }
      
      //
      /// describes the strong form of the PDE
      //
      virtual void L(DoubleVector& b, const FemFunction& U) const=0;
      
      //
      /// describes the stabilization term of the PDE;
      /// can be chosen as -L^
      //
      virtual void S(DoubleMatrix& A, const FemFunction& U, const TestFunction& N) const=0;

      //
      /// describes the strong derivative of the PDE
      //
      virtual void LMatrix(DoubleMatrix& A, const FemFunction& U, const TestFunction& M) const=0;

      //
      /// describes the derivative of the stabilization term S;
      //
      virtual void SMatrix(DoubleVector& b, const FemFunction& U, const FemFunction& M, const FemFunction& N) const {};
  };
}

/*-----------------------------------------*/

#endif
