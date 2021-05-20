
#include "boundaryfunctional.h"
#include "dirichletdatabycolor.h"
#include "local.h"
#include "loop.h"
#include "residualfunctional.h"
#include "weightedpointfunctional.h"

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
    __cols.insert(84);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
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
    __cols.insert(84);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};

/*---------------------------------------------------*/

int
main(int argc, char** argv)
{
  ParamFile pf("box3d.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor3d Problem3d;
  Problem3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("fsi", &Problem3d);

  FunctionalContainer FC3d;
  WeightedPointFunctional Ux;
  WeightedPointFunctional Uy;
  WeightedPointFunctional Uz;
  vector<Vertex3d> v1;
  v1.push_back(Vertex3d(0.45, 0.15, 0.15));

  vector<int> cx;
  cx.push_back(4);
  vector<int> cy;
  cy.push_back(5);
  vector<int> cz;
  cz.push_back(6);

  vector<double> weigh;
  weigh.push_back(1.0);
  /*Ux.BasicInit(v1,cx,weigh);
  Uy.BasicInit(v1,cy,weigh);
  Uz.BasicInit(v1,cz,weigh);

  Drag drag;
  Lift lift;

  FC3d.AddFunctional("ux",&Ux);
  FC3d.AddFunctional("uy",&Uy);
  FC3d.AddFunctional("uz",&Uz);
  FC3d.AddFunctional("drag",&drag);
  FC3d.AddFunctional("lift",&lift);
  */

  Loop<3> loop;

  loop.BasicInit(&pf, &PC3d, &FC3d);
  loop.run("fsi");

  return 0;
}
