#ifndef  __NavierStokes3d_h
#define  __NavierStokes3d_h

#include  "navierstokes2d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class NavierStokes3d : public NavierStokes2d
{
  protected:
  
  double Laplace(const TestFunction& U, 
		 const TestFunction& N) const
    {
      return U.x()*N.x() + U.y()*N.y() + U.z()*N.z();
    }
  
  double Convection(const std::vector<TestFunction>& U, 
		    const TestFunction& N) const
    {
      return U[1].m()*N.x() + U[2].m()*N.y() + U[3].m()*N.z();
    }
  double Divergence(const std::vector<TestFunction>& U) const
    {
      return U[1].x() + U[2].y() + U[3].z();
    }
  
  public:

  ~NavierStokes3d();
  NavierStokes3d();
  NavierStokes3d(const ParamFile* pf);

  std::string GetName() const;

  int    GetNcomp  () const { return 4; }

  DoubleVector GetViscosities() const;

  void SetTimePattern(TimePattern& P) const;

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  void point(double h, const FemFunction& U, const Vertex3d& v) const 
    { 
      _h = h;
    }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
  void MatrixLoop(EntryMatrix& A, const FemFunction& U, const FemFunction& M, const FemFunction& N) const;
};
}

#endif
