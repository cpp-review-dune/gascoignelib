#ifndef  __NavierStokes3d_h
#define  __NavierStokes3d_h

#include  "navierstokes.h"

/*-----------------------------------------*/

class NavierStokes3d : public NavierStokes
{
  protected:
  
  double Laplace(const Gascoigne::TestFunction& U, 
		 const Gascoigne::TestFunction& N) const
    {
      return U.x()*N.x() + U.y()*N.y() + U.z()*N.z();
    }
  
  double Convection(const std::vector<Gascoigne::TestFunction>& U, 
		    const Gascoigne::TestFunction& N) const
    {
      return U[1].m()*N.x() + U[2].m()*N.y() + U[3].m()*N.z();
    }
  double Divergence(const std::vector<Gascoigne::TestFunction>& U) const
    {
      return U[1].x() + U[2].y() + U[3].z();
    }
  
  public:

  ~NavierStokes3d();
  NavierStokes3d();
  NavierStokes3d(const Gascoigne::ParamFile* pf);

  std::string GetName() const;

  int    ncomp  () const { return 4; }

  nvector<double> GetViscosities() const;

  void SetTimePattern(TimePattern& P) const;

  void OperatorStrong(Vector& b, const Gascoigne::FemFunction& U) const;

  void point(double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const 
    { 
      _h = h;
    }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;

  void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;
};

#endif
