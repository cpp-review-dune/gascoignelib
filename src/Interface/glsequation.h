#ifndef  __GlsEquation_h
#define  __GlsEquation_h

#include "equation.h"

//////////////////////////////////////////////
//
///@brief
/// Interface class for  Elements

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

class GlsEquation : public virtual Equation
{
public:

  GlsEquation() {};
  ~GlsEquation() {};

  //
  /// computation of stabilization parameters
  //
  virtual void glspoint
    (double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const { assert(0);}
  //
  /// computation of stabilization parameters
  //
  virtual void glspoint
    (double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const { assert(0);}
 
  virtual void glspoint(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Vertex2d& v) const {
    glspoint(h,U,v);
  }
  virtual void glspoint(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Vertex3d& v) const {
    glspoint(h,U,v);
}
  //
  /// computation of stabilization parameters for the matrix
  //
  virtual void glspointmatrix(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Vertex2d& v) const {
    glspoint(h,U,Q,v);
  }
  virtual void glspointmatrix(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Vertex3d& v) const {
    glspoint(h,U,Q,v);
  }
  //
  /// describes the strong form of the PDE
  //
  virtual void L(nvector<double>& dst, const Gascoigne::FemFunction& U) const=0;
  //
  /// describes the stabilization term of the PDE;
  /// can be chosen as -L^
  //
  virtual void S(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const=0;

  //
  /// describes the strong derivative of the PDE
  //
  virtual void LMatrix(nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const=0;

  //
  /// describes the derivative of the stabilization term S;
  //
  virtual void SMatrix(nvector<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::FemFunction& M, const Gascoigne::FemFunction& N) const {}
};

/*-----------------------------------------*/

#endif
