#include  "problemdescriptor1.h"
#include  "stdloop.h"
#include  "starter.h"
#include  "meshagent.h"
#include  "boundaryfunction.h"

/* ----------------------------------------- */
class RunderKreis : public BoundaryFunction<2>
{
  double squareradius;
  Vertex2d center;
public :
  void BasicInit(Vertex2d c, double r) {
    center = c; 
    squareradius = r;
  }

  double operator()(const Vertex2d& c) const {
    double r = - squareradius;
    for (int i=0; i<2; i++)
      {
	double dx = c[i]-center[i];
	r += dx * dx;
      }
    return r;
  }
};

RunderKreis RK;

/* ----------------------------------------- */
class LocalLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile* paramfile) {
  GetMeshAgentPointer() = new MeshAgent;

  double r = 0.25;
  Vertex2d v(2.,2.);
  RK.BasicInit(v,r);

  int prerefine=1;
  string inpname("nsbench4.inp");
  std::map<int,BoundaryFunction<2>* > shapes;
  shapes[80] = &RK;
  dynamic_cast<MeshAgent*>(GetMeshAgent())->BasicInit(inpname, prerefine,shapes);

  StdLoop::BasicInit(paramfile);
  }
};

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  Starter S(argc, argv, "bench.param");

  ProblemDescriptor1 LPD;
  LPD.BasicInit(S.GetParamFile());

  LocalLoop loop;
  loop.BasicInit(S.GetParamFile());
  loop.run(&LPD);

  return 0;
}
