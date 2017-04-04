#include  "local.h"
#include  "aleloop.h"
#include  "weightedpointfunctional.h"


using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/


#include  "domainfunctional.h"
#include  "residualfunctional.h" 
#include  "boundaryfunctional.h" 
#include  "dirichletdatabycolor.h"


bool __ADJOINT;  
bool __ADJOINT1;  
double __EXACT//  = 0.00943419204958;
= 2.27e-5;

bool __ESTIMATE = false;



int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor3d LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("primal", &LPD);

  FunctionalContainer FC;


  /////// forces
  LocalDomainFunctionals_FlagForce drag("drag-3d");
  drag.ExactValue() = 1.327;
  DirichletDataByColor    drag_dd(drag.GetComps(), drag.GetColors(), drag.GetScales());
  DualProblemDescriptor3d drag_dual(&paramfile,"drag-3d", NULL, &drag_dd);  
  FC.AddFunctional("drag-3d",&drag);
  PC.AddProblem("drag-3d",&drag_dual);

  // LocalDomainFunctionals_FlagForce lift("lift-3d");
  // lift.ExactValue() = 0.0;
  // DirichletDataByColor    lift_dd(drag.GetComps(), drag.GetColors(), drag.GetScales());
  // DualProblemDescriptor3d lift_dual(&paramfile,"lift-3d", NULL, &drag_dd);  
  // FC.AddFunctional("lift-3d",&lift);
  // PC.AddProblem("lift-3d",&lift_dual);


  Vertex3d v(0.45,0.15,0.15);
  vector<Vertex3d> vv; vv.push_back(v);
  vector<int> cps(1); 
  vector<double> w(1,1.0);
  ZeroDirichletData ZeroDD;

  WeightedPointFunctional wx; cps[0] = 4;
  wx.BasicInit(vv,cps,w);
  wx.ExactValue() = 5.924e-5;
  
  WeightedDiracRightHandSide wx_rhs; wx_rhs.BasicInit(&wx);
  DualProblemDescriptor3d wx_dual(&paramfile,"wx", &wx_rhs, &ZeroDD);
  FC.AddFunctional("wx", &wx);
  PC.AddProblem("wx", &wx_dual);

  // WeightedPointFunctional wy; cps[0] = 5;
  // wy.BasicInit(vv,cps,w);
  // wy.ExactValue() = 0.0;
  // WeightedDiracRightHandSide wy_rhs; wy_rhs.BasicInit(&wy);
  // DualProblemDescriptor3d wy_dual(&paramfile,"wy", &wy_rhs, &ZeroDD);
  // FC.AddFunctional("wy", &wy);
  // PC.AddProblem("wy", &wy_dual);

  // WeightedPointFunctional wz; cps[0] = 6;
  // wz.BasicInit(vv,cps,w);
  // wz.ExactValue() = 0.0;
  // WeightedDiracRightHandSide wz_rhs; wz_rhs.BasicInit(&wz);
  // DualProblemDescriptor3d wz_dual(&paramfile,"wz", &wz_rhs, &ZeroDD);
  // FC.AddFunctional("wz", &wz);
  // PC.AddProblem("wz", &wz_dual);
  


  
  AleLoop loop;
  
  loop.BasicInit(&paramfile, &PC, &FC);
  loop.run("primal");

  return 0;
}
