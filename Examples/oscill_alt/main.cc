#include <omp.h>
#include "local.h"
#include "loop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include "residualfunctional.h"
#include "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;

class Drag : public virtual ResidualFunctional {
  std::string GetName() const {
    return "drag";
  }

public:
  Drag() {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(80);
    __DD= new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};
class Lift : public virtual ResidualFunctional {
  std::string GetName() const {
    return "lift";
  }

public:
  Lift() {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(80);
    __DD= new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};

/*---------------------------------------------------*/

int main(int argc, char** argv) {
  auto start= omp_get_wtime();
#ifdef BDF
  ParamFile pf("fsi-3.param");
#else
  ParamFile pf("fsi-wm.param");
#endif
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor2d Problem2d;
  Problem2d.BasicInit(&pf);

  DivergenceDescriptor2d Divergence2d;
  Divergence2d.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("fsi", &Problem2d);
  PC2d.AddProblem("div", &Divergence2d);

  FunctionalContainer     FC2d;
  WeightedPointFunctional Px, Py;
  vector<Vertex2d>        v2d;
  v2d.push_back(Vertex2d(0.6, 0.2));
  vector<int> cx, cy;
  cx.push_back(3);
  cy.push_back(4);
  vector<double> weigh;
  weigh.push_back(1.0);
  Px.BasicInit(v2d, cx, weigh);
  Py.BasicInit(v2d, cy, weigh);
  Drag drag;
  Lift lift;

  FC2d.AddFunctional("ux", &Px);
  FC2d.AddFunctional("uy", &Py);
  FC2d.AddFunctional("drag", &drag);
  FC2d.AddFunctional("lift", &lift);

  Loop<2> loop;

  loop.BasicInit(&pf, &PC2d, &FC2d);
  loop.run("fsi");
  auto end        = omp_get_wtime();
  auto serial_time= end - start;
  std::cerr << "\nElapsed time is\t" << serial_time << ".\n";
  return 0;
}
