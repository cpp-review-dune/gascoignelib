#ifndef  __RightHandSideByEquation_h
#define  __RightHandSideByEquation_h

#include  "domainrighthandside.h"
#include  "exactsolution.h"
#include  "equation.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class RightHandSideByEquation : public DomainRightHandSide
{
protected:

  const Equation*      _EQ;
  const ExactSolution* _ES;

public:

  RightHandSideByEquation(const Equation* eq, const ExactSolution* es)
    : DomainRightHandSide(), _EQ(eq), _ES(es) { assert(es); assert(eq);}
  
  std::string GetName() const { return "RightHandSideByEquation";} 
  int GetNcomp() const { return _EQ->GetNcomp();}

  double operator()(int c, const Vertex2d& v)const 
    {
      int n = _EQ->GetNcomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].m() = (*_ES)(i,v);
	  U[i].x() = _ES->x(i,v);
	  U[i].y() = _ES->y(i,v);
	  U[i].D() = _ES->xx(i,v)+_ES->yy(i,v);
	}
      _EQ->OperatorStrong(b,U);
      if (GetTimeStep()>0.)
	{
	  double eps = 1.e-6;
	  DoubleVector ut(n,0.);
	  double time = GetTime();
	 _ES->SetTime(time+0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time,GetTimeStep());
	  b.add(1./eps,ut);
	}
      return b[c];
    }
  double operator()(int c, const Vertex3d& v)const 
    {
      int n = _EQ->GetNcomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].m() = (*_ES)(i,v);
	  U[i].x() = _ES->x(i,v);
	  U[i].y() = _ES->y(i,v);
	  U[i].z() = _ES->z(i,v);
	  U[i].D() = _ES->xx(i,v)+_ES->yy(i,v)+_ES->zz(i,v);
	}
      _EQ->OperatorStrong(b,U);
      if (GetTimeStep()>0.)
	{
	  double eps = 1.e-6;
	  DoubleVector ut(n,0.);
	  double time = GetTime();
	 _ES->SetTime(time+0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps,GetTimeStep());
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time,GetTimeStep());
	  b.add(1./eps,ut);
	}
      return b[c];
    }

};
}

#endif
