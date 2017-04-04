#include  "problem.h"
#include  "loop.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  if (argc<4)
    {
      cerr << "Theta test.param  NITER  THETA " << endl;
      abort();
    }
  int  niter = atoi(argv[2]);
  double theta = atof(argv[3]);



  //////////////
  // Functional
  /////////////
  WeightedPointFunctional P;
  vector<Vertex2d> vv;vv.push_back(Vertex2d(0.0,-0.5));
  vector<double>   ww;ww.push_back(1.0);
  vector<int>      cc;cc.push_back(0);
  P.BasicInit(vv,cc,ww);
  
  MyFunc F1;
  
  FunctionalContainer FC;
  FC.AddFunctional("end", &F1);
  FC.AddFunctional("time", &F1);

  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("test", &LPD);
  

  
  /////////////
  // Loop
  /////////////
  ThetaLoop loop;
  loop.BasicInit(&paramfile, &PC, &FC);
  
  loop.run_exact("test", niter, theta);
  //  loop.run("test", niter, theta);

  return 0;
}
