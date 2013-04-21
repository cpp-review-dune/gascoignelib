#include  "boundaryfunction.h"
#include <cmath>

using namespace std;

/*---------------------------------------------------*/

namespace Tsuchimikado
{
  template<int DIM>
  void BoundaryFunction<DIM>::grad(Gascoigne::Vertex<DIM>& dst, const Gascoigne::Vertex<DIM>& src) const 
  {
    double eps = 1e-6;
  
    for(int i=0;i<DIM;i++)
      {
	Gascoigne::Vertex<DIM> cl(src), cr(src);
	cl[i] -= eps;
	cr[i] += eps;
	dst[i] = ((*this)(cr)-(*this)(cl))/(2.*eps);
      }
  }

  /*---------------------------------------------------*/

  template<int DIM>
  void BoundaryFunction<DIM>::newton(Gascoigne::Vertex<DIM>& dst) const
  {
    int    maxi = 10;
    double tol  = 1.e-12;

    Gascoigne::Vertex<DIM> z;

    grad(z,dst);
    double res = (*this)(dst);

    for (int i=0; (i<maxi) && (fabs(res)>tol) ; i++)
      {
	Gascoigne::Vertex<DIM> zz;
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

  template class BoundaryFunction<2>;
  template class BoundaryFunction<3>;
}
