#ifndef  __RightHandSideDataByEquation_h
#define  __RightHandSideDataByEquation_h

#include  "righthandsidedata.h"
#include  "exactsolution.h"
#include  "equation.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class RightHandSideDataByEquation : public RightHandSideData
{
protected:

  const Equation*      _EQ;
  const ExactSolution* _ES;

public:

  RightHandSideDataByEquation(const Equation* eq, const ExactSolution* es)
    : RightHandSideData(), _EQ(eq), _ES(es) { assert(es); assert(eq);}
  
  std::string GetName() const { return "laplace";} 
  int GetNcomp() const { return _EQ->ncomp();}

  double operator()(int c, const Vertex2d& v)const 
    {
      int n = _EQ->ncomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].n() = 0.;
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
	 _ES->SetTime(time+0.5*eps);
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps);
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time);
	  b.add(1./eps,ut);
	}
      return b[c];
    }
  double operator()(int c, const Vertex3d& v)const 
    {
      int n = _EQ->ncomp();
      DoubleVector b(n,0.);
      FemFunction U(n);
      for (int i=0; i<n; i++)
	{
	  U[i].n() = 0.;
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
	 _ES->SetTime(time+0.5*eps);
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] = (*_ES)(i,v);
	    }
	  _ES->SetTime(time-0.5*eps);
	  for (int i=0; i<n; i++)
	    {	  
	      ut[i] -= (*_ES)(i,v);
	    }
	  _ES->SetTime(time);
	  b.add(1./eps,ut);
	}
      return b[c];
    }

};
}

#endif
