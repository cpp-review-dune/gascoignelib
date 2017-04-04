
#include  "local.h"
#include  "loop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include "functionals.h"


using namespace Gascoigne;
using namespace std;

extern nvector<double> exakt;





/*---------------------------------------------------*/



int main(int argc, char** argv)
{
  ParamFile pf("bench.param");
  assert(argc==1);

  // FUNCTIONAL
  //  exakt.push_back(16.97959184); // Mittelwert
  //
  //  exakt.push_back(20.565374463129671767); 
  //  exakt.push_back(0.0103708331251164982);
  //  exakt.push_back(25.1595); // v end
  //  exakt.push_back(0.056507); // lift end


  //  exakt.push_back(-0.023368109455691368076);

  
  ///////// Beispiel 2: (Testcase 1)
  exakt.push_back(2702.11801235531766);// +/-
  //  exakt.push_back(-4.4711914866058553031);//
  




  FunctionalContainer FC2d;
  VMEAN vmean;
  FC2d.AddFunctional("0 vmean", &vmean);
  //  Blift lift(&pf);
  //  FC2d.AddFunctional("0 lift", &lift);
  //LiftBRHS lift_rhs(&pf);
  
  //  Pdrop pdrop;
  

  // PROBLEM
  
  ProblemDescriptor2d primal;
  primal.BasicInit(&pf);
  AdjointProblemDescriptor2d adjoint;
  VMEANRhs vmean_rhs;

  adjoint.BasicInit(&pf,vmean_rhs);
  //adjoint.BasicInit(&pf,pdrop);
  //  adjoint.BasicInit(&pf,lift_rhs);


  ProblemContainer PC2d;
  PC2d.AddProblem("primal", &primal);
  PC2d.AddProblem("adjoint", &adjoint);

  Loop loop;
  loop.BasicInit(&pf,&PC2d,&FC2d);
  
  //Berechnet die Fehlerschatzer verfeinert aber trotzdem global
  //loop.run("primal");
  
  //Berechnet die Fehlerschatzer und verfeinert entsprechend
  loop.run_fst_adaptive("primal");

  return 0;
}
