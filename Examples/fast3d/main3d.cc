
#include  "local.h"
#include  "loop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;


class Drag : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "drag"; } 
public:
  Drag()
  {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(80);   
    __cols.insert(81);   
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());

  }
};
class Lift : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "drag"; } 
public:
  Lift()
  {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(80);   
    __cols.insert(81);   
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());

  }
};




/*---------------------------------------------------*/



int main(int argc, char** argv)
{
  ParamFile pf("fsi-3d.param");
  if (argc==2)
    pf.SetName(argv[1]);

  
  ProblemDescriptor3d Problem3d;
  Problem3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("fsi", &Problem3d);

  FunctionalContainer FC3d;
  WeightedPointFunctional U1x,U1y,U1z;
  WeightedPointFunctional U2x,U2y,U2z;
  WeightedPointFunctional U3x,U3y,U3z;
  vector<Vertex3d> v1;  v1.push_back(Vertex3d(0.90,0.20,0.10));
  vector<Vertex3d> v2;  v2.push_back(Vertex3d(0.90,0.20,0.20));
  vector<Vertex3d> v3;  v3.push_back(Vertex3d(0.90,0.20,0.30));

  vector<int> cx; cx.push_back(4);
  vector<int> cy; cy.push_back(5);
  vector<int> cz; cz.push_back(6);

  vector<double> weigh; weigh.push_back(1.0);
  U1x.BasicInit(v1,cx,weigh); 
  U1y.BasicInit(v1,cy,weigh); 
  U1z.BasicInit(v1,cz,weigh); 

  U2x.BasicInit(v2,cx,weigh); 
  U2y.BasicInit(v2,cy,weigh); 
  U2z.BasicInit(v2,cz,weigh); 

  U3x.BasicInit(v3,cx,weigh); 
  U3y.BasicInit(v3,cy,weigh); 
  U3z.BasicInit(v3,cz,weigh); 

  Drag drag;
  Lift lift;
  
  FC3d.AddFunctional("u1x",&U1x);
  FC3d.AddFunctional("u1y",&U1y);
  FC3d.AddFunctional("u1z",&U1z);

  FC3d.AddFunctional("u2x",&U2x);
  FC3d.AddFunctional("u2y",&U2y);
  FC3d.AddFunctional("u2z",&U2z);

  FC3d.AddFunctional("u3x",&U3x);
  FC3d.AddFunctional("u3y",&U3y);
  FC3d.AddFunctional("u3z",&U3z);
  FC3d.AddFunctional("drag",&drag);
  FC3d.AddFunctional("lift",&lift);
  
  Loop<3> loop;

  loop.BasicInit(&pf,&PC3d,&FC3d);
  loop.run("fsi");

  return 0;
}
