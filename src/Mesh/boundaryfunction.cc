#include  "boundaryfunction.h"

using namespace std;

/*---------------------------------------------------*/

template<int DIM>
void BoundaryFunction<DIM>::grad
(numfixarray<DIM,double>& dst, const numfixarray<DIM,double>& src) const 
{
  double eps = 1e-6;
  
  for(int i=0;i<DIM;i++)
    {
      numfixarray<DIM,double> cl(src), cr(src);
      cl[i] -= eps;
      cr[i] += eps;
      dst[i] = ((*this)(cr)-(*this)(cl))/(2.*eps);
    }
}

/*---------------------------------------------------*/

template<int DIM>
void BoundaryFunction<DIM>::newton(Vector& dst) const
{
  int    maxi = 10;
  double tol  = 1.e-12;

  Vector z;

  grad(z,dst);
  double res = (*this)(dst);

  for (int i=0; (i<maxi) && (fabs(res)>tol) ; i++)
    {
      Vector zz;
      grad(zz,dst);
      double bgrad = z*zz;

      if (fabs(bgrad)<=1.e-15)
	{
	  cerr << "BoundaryFunction<DIM>::newton()\n";
	  cerr << "Grad=0 in boundary_newton (res= "<< res << " )\n";
	  cerr << "iter " << i << endl;
	  abort();
	}
      res /= bgrad;
      dst.add(-res,z);
      res = (*this)(dst);
    }
  if (fabs(res)>tol)
    {
      cerr << "BoundaryFunction<DIM>::newton()\n";
      cerr << "No Convergence in boundary_newton (res= "<< res << " )\n";
      abort();
    }
}

/*---------------------------------------------------*/

template BoundaryFunction<2>;
template BoundaryFunction<3>;
