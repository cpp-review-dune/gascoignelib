#include "meshagent.h"
#include "stdsolver.h"
#include "stdmultilevelsolver.h"
#include "gascoigne.h"
#include "paramfile.h"
#include "q12d.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

class LocalSolver : public StdSolver
{
public:

  LocalSolver() : StdSolver() {};
  ~LocalSolver() {}

  double Integral(const GlobalVector& u) const
    {
      HNAverage(u);
      assert(GetDiscretization()->HNZeroCheck(u)==0);
      
      nvector<double> dst = _PF.IntegrateVector(u);
      HNZero(u);
      return dst[0];
    }
  DiscretizationInterface* NewDiscretization(int dimension, const string& discname)
    {
      return new Q12d;
    }
  void NewMesh(int l, const MeshInterface* MP)
    {
      StdSolver::NewMesh(l,MP);
      GetDiscretization()->InitFilter(_PF);
    }
  void BasicInit(int level, const ParamFile* paramfile, const MeshInterface* MP)
    {
      StdSolver::BasicInit(level,paramfile,MP);
    }
};

/*---------------------------------------------------*/

double f(double x, double t)
{
  double xx = 0.5 + 0.25*cos(2*M_PI*t);
  double yy = 0.5 + 0.25*sin(2*M_PI*t);
  double u = 1./(1.+(x-xx)*(x-xx)+yy*yy);
  double val = u*u*yy;
  return val;
}
/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }

  int niter = 10;

  MeshAgent M;
  M.BasicInit(&paramfile);

  GlobalVector u(1);

  for (int iter=1; iter<=niter; iter++)
    {
      M.global_refine(1);

      LocalSolver S;

      const MeshInterface* MI = M.GetMesh(0);
      S.BasicInit(iter,&paramfile,MI);
      S.NewMesh(iter,MI);
      
      u.resize(MI->nnodes());
      
      for (int i=0; i<u.n(); i++)
	{
	  double x = M.GetMesh()->vertex2d(i).x();
	  double y = M.GetMesh()->vertex2d(i).y();
	  u(i,0) = f(x,y);
	}
      
      double val = S.Integral(u);
      cout.precision(16);
      cout << iter << "\t" << u.n() << "\t" << val << endl;
    }
}
