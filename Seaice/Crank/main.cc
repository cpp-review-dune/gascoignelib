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
/*

class Temp  : public virtual BoundaryFunctional
{

  
public:


  Temp(const ParamFile* pf)
  {
    
  }

  std::string GetName() const {return "Temp";}

  double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
  {
    double c= 0.249702115724673;
    
    if(v.y()==0)
      return ( 1/c * U[0].y() );
    
    return 0;
  }
};


*/

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
   double x=v.x()-0.25;
  double  y=v.y()-0.25;
   

     return  exp(-10*(x*x+y*y))*U[0].m();

    
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

  DualHelpProblem DHP;
  DHP.BasicInit(&paramfile);

  // NonlinearProblem NON;
  // NON.BasicInit(&paramfile);


  ProblemContainer PC;
  PC.AddProblem("LaplaceT", &LPD);
  PC.AddProblem("dp",    &DP);
  PC.AddProblem("dhp",&DHP);
  //  PC.AddProblem("Non",&NON);
  
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
