#include <omp.h>
#include "local.h"
#include "loop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include "residualfunctional.h"
#include "resfunctional.h"
#include "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv) {
  auto      start= omp_get_wtime();
  ParamFile pf("fsi_box3d.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor3d Problem3d;
  Problem3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("fsi", &Problem3d);

  FunctionalContainer     FC3d;
  WeightedPointFunctional Ux;
  WeightedPointFunctional Uy;
  WeightedPointFunctional Uz;
  vector<Vertex3d>        v1;
  v1.push_back(Vertex3d(0.45, 0.15, 0.15));

  vector<int> cx;
  cx.push_back(4);
  vector<int> cy;
  cy.push_back(5);
  vector<int> cz;
  cz.push_back(6);

  vector<double> weigh;
  weigh.push_back(1.0);
  Ux.BasicInit(v1, cx, weigh);
  Uy.BasicInit(v1, cy, weigh);
  Uz.BasicInit(v1, cz, weigh);

  Drag3d drag;
  Lift3d lift;

  FC3d.AddFunctional("ux", &Ux);
  FC3d.AddFunctional("uy", &Uy);
  FC3d.AddFunctional("uz", &Uz);
  FC3d.AddFunctional("drag", &drag);
  FC3d.AddFunctional("lift", &lift);

  Loop<3> loop;

  loop.BasicInit(&pf, &PC3d, &FC3d);
  loop.run("fsi");
  auto end        = omp_get_wtime();
  auto serial_time= end - start;
  std::cerr << "\nElapsed time is\t" << serial_time << ".\n";
  return 0;
}
