#include  "problem.h"
#include  "loop.h"
#include "boundaryfunctional.h"
#include  "filescanner.h"
#include  "stdloop.h"
#include "gascoignemesh2d.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "simplematrix.h"
using namespace Gascoigne;
using namespace std;

extern double TIME,DT;

class Temp  : public virtual DomainFunctional
{
  
public:
 

  
  Temp(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
  }
  std::string GetName() const {return "Ellipse";}
  double J(const FemFunction& U, const Vertex2d& v) const
  { 
    if(TIME>3.0 && TIME<3.5)
    { double x=v.x()-0.5; double y=v.y()-0.5;
      return  (U[0].x()-U[1].y())*exp(-x*x-y*y);}
    
    return 0;
  }
};






int main(int argc, char** argv)
{
  // parameter file as argument
  ParamFile paramfile("gascoigne.param");
  if(argc>=2) 
    paramfile.SetName(argv[1]);
  

  // define problem
  LaplaceTProblem LPD;
  LPD.BasicInit(&paramfile);

    DualProblem DP;
  DP.BasicInit(&paramfile);

 


  ProblemContainer PC;
  PC.AddProblem("LaplaceT", &LPD);
  PC.AddProblem("dp",    &DP);
 
 
  
  FunctionalContainer FC;
  Temp T(&paramfile);


  FC.AddFunctional("0 Temp", &T);
  
  // loop for program control
  Loop loop;
  loop.BasicInit(&paramfile, &PC, &FC);

  // start solving
  loop.run("LaplaceT");

  return 0;
}