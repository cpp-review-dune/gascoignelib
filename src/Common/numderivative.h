#ifndef __numderivative_h
#define __numderivative_h

template <class C, class MAT, class VEC>
void numderivative(C& application, MAT& M, const VEC& x, double eps=1.e-4)
{
  int m = x.size();

  VEC up(x), u(x), xp(x);
  up.zero(); u.zero();

//  double ieps = 1./eps;
  
  application.f(u,x);
  
  for (int i=0; i<m ; i++)
    {
      //xp[i] += eps;
      xp[i] *= 1.+eps;
      application.f(up,xp);
      for (int j=0; j<m ; j++)
	{
	  //M(j,i) = (up[j]-u[j])*ieps;
	  M(j,i) = (up[j]-u[j])/(xp[i]*eps);
	}
      //xp[i] -= eps;
      xp[i] /= 1.+eps;
    }
}

#endif
