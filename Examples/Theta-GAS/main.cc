#include  "problem.h"
#include  "loop.h"
#include  "weightedpointfunctional.h"


using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

extern double exact_time;

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
  



  // loop.run_exact("test", 20, 0.5);
  // //loop.run_adaptive("test", niter, theta);
  

  

  
  

  //   exact_time = 0.256743392; // p=4
  //   exact_time =0.265484037 ; // p=2
  //   exact_time = 0.2665432686; // p=1.5
  
  //   exact_time = 0.3712646796; //  (p=1.5, extra beta)
  //               0.3712615584
  //   exact_time = 0.2355601141;//0.2355601393;
     
  //   exact_time = 14.48498326;
  
  //Example 1
    exact_time = -2.3632006242; //p=2 ,eps=0.01
  //Example 2   
  //   exact_time = 0.256743392; // p=4 ,eps=0.001
  //   exact_time =0.265484037 ; // p=2 ,eps=0.001
  //   exact_time = 0.2665432686; // p=1.5,eps=0.001
  
  
     
     
     
     
     
     
     
     
     
     int iter=0;
     
     for (int M=24;M<1000;M*=2)
     {
       cout<<" Start run theta with M= "<<M<<"and  theta="<<theta<<endl;
       loop.SetReflevel(iter);
       loop.run_theta("test", M, theta);
       iter++;
     }

     //   loop.run("test", M, theta);
     
    //loop.run_adaptive_theta("test", niter, theta);
     //loop.run_adaptive("test", niter, theta);

 
  return 0;
}
