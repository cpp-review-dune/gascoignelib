/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "boundaryfunction.h"
#include "gascoignehash.h"
#include "gascoignemesh2d.h"
#include "local.h"
#include "meshagent.h"
#include "solvers.h"
#include "stdloop.h"
#include "usefullfunctionsbd.h"
#include "vertex.h"

namespace Gascoigne {

class Cyl : public BoundaryFunction<3>
{
  double squareradius;
  Vertex2d center;

public:
  std::string GetName() const { return "Cyl"; }

  void BasicInit(Vertex2d c, double r)
  {
    center = c;
    squareradius = r;
  }

  double operator()(const Vertex3d& c) const
  {
    double r = -squareradius;
    for (int i = 0; i < 2; i++) {
      double dx = c[i] - center[i];
      r += dx * dx;
    }
    return r;
  }
};

class Cir : public BoundaryFunction<2>
{
  double squareradius;
  Vertex2d center;

public:
  std::string GetName() const { return "Cir"; }

  void BasicInit(Vertex2d c, double r)
  {
    center = c;
    squareradius = r;
  }

  double operator()(const Vertex2d& c) const
  {
    double r = -squareradius;
    for (int i = 0; i < 2; i++) {
      double dx = c[i] - center[i];
      r += dx * dx;
    }
    return r;
  }
};

/*---------------------------------------------------*/

class MA3d : public MeshAgent
{
protected:
  Cyl RK;

public:
  MA3d()
    : MeshAgent()
  {
    double r2 = 0.012 * 0.012;
    Vertex2d v(0.0, 0.0);
    RK.BasicInit(v, r2);
    AddShape(4, &RK);
  }
};

class MA2d : public MeshAgent
{
protected:
  Cir RK;

public:
  MA2d()
    : MeshAgent()
  {
    double r2 = 0.05 * 0.05;
    Vertex2d v(0.2, 0.2);
    RK.BasicInit(v, r2);
    AddShape(80, &RK);
    AddShape(81, &RK);
  }
};

template<int DIM>
class Loop : public StdLoop
{

public:
  void BasicInit(const ParamFile* paramfile,
                 const ProblemContainer* PC,
                 const FunctionalContainer* FC)
  {
    if (DIM == 2)
      GetMeshAgentPointer() = new MA2d;
    if (DIM == 3)
      GetMeshAgentPointer() = new ProjectionOnFineMeshAgent;
    GetMultiLevelSolverPointer() = new FSIMultiLevelSolver<DIM>;

    StdLoop::BasicInit(paramfile, PC, FC);
  }

  void run(const std::string& problemlabel);
};

} // namespace Gascoigne

/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
