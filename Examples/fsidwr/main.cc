#include  "local.h"
#include  "aleloop.h"
#include  "weightedpointfunctional.h"
#include  "chi.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/



#include  "domainfunctional.h"
#include  "residualfunctional.h" 
#include  "boundaryfunctional.h" 
#include  "dirichletdatabycolor.h"




bool __ADJOINT;  
bool __ADJOINT1;  
double __EXACT= 2.27e-5;//  = 0.00943419204958;





bool __ESTIMATE;

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh1.param");
  if(argc>=2) {
    paramfile.SetName(argv[1]);
  }
  
  /////////////
  // Equation
  /////////////
  ProblemDescriptor LPD;
  LPD.BasicInit(&paramfile);

  ProblemContainer PC;
  PC.AddProblem("primal", &LPD);

  

  //////////////////////////////////////////////////
  ////////////////////////////////////////////////// Functionals and Adjoint Problems
  //////////////////////////////////////////////////
  ZeroDirichletData ZeroDD;

  FunctionalContainer FC;

  ///// Total derivation 'df'
  DF df(&paramfile);
  df.ExactValue() =0.0004493943165;

  // FC.AddFunctional("deri", &df);
  // DF_RHS<2> df_rhs(&paramfile);
  // DualProblemDescriptor  df_dual(&paramfile, "deri", &df_rhs, &ZeroDD);
  // PC.AddProblem("deri", &df_dual);
  
  ///// Point evaluations 'wx', 'wy'


  Vertex2d v(0.6,0.2);
  vector<Vertex2d> vv; vv.push_back(v);
  vector<int> cps(1); 
  vector<double> w(1,1.0);

  // WeightedPointFunctional wx; cps[0] = 3;
  // wx.BasicInit(vv,cps,w);
  // wx.ExactValue() = 2.2685e-5;
  // FC.AddFunctional("wx", &wx);
  // WeightedDiracRightHandSide wx_rhs; wx_rhs.BasicInit(&wx);
  // DualProblemDescriptor wx_dual(&paramfile,"wx", &wx_rhs, &ZeroDD);
  // PC.AddProblem("wx", &wx_dual);
  // assert(wx_rhs.GetPoints2d().size()==1);
  
  // WeightedPointFunctional wy; cps[0] = 4;  
  // wy.BasicInit(vv,cps,w);
  // wy.ExactValue() = 8.186209e-4;
  // FC.AddFunctional("wy", &wy);
  // WeightedDiracRightHandSide wy_rhs; wy_rhs.BasicInit(&wy);  
  // DualProblemDescriptor wy_dual(&paramfile,"wx", &wy_rhs, &ZeroDD);
  // PC.AddProblem("wy", &wy_dual);
  // assert(wy_rhs.GetPoints2d().size()==1);



  /////// forces
  LocalDomainFunctionals_FlagForce drag("drag");
  drag.ExactValue() = 14.29395;

  ///// nur fluid!!!
  drag.ExactValue() = 0.0;
  
    FC.AddFunctional("drag", &drag);
  DirichletDataByColor drag_dd(drag.GetComps(), drag.GetColors(), drag.GetScales());
  DualProblemDescriptor drag_dual(&paramfile,"drag", NULL, &drag_dd);  
  PC.AddProblem("drag", &drag_dual);
  


  // LocalDomainFunctionals_FlagForce lift("lift");
  // lift.ExactValue() = 0.764761;
  // FC.AddFunctional("lift", &lift);
  // DirichletDataByColor lift_dd(lift.GetComps(), lift.GetColors(), lift.GetScales());
  // DualProblemDescriptor lift_dual(&paramfile,"lift", NULL, &lift_dd);
  // PC.AddProblem("lift", &lift_dual);

  

  
  AleLoop loop;
  
  loop.BasicInit(&paramfile, &PC, &FC);
  loop.run("primal");

  return 0;
}
