#ifndef __BenchMarkMeshAgent_h
#define __BenchMarkMeshAgent_h

#include "runderkreis.h"
#include "meshagent.h"

/*---------------------------------------------------*/

class BenchMarkMeshAgent : public MeshAgent
{
 protected:
  
  RunderKreis RK;

 public:
  
  BenchMarkMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      int prerefine=1;
      string inpname("nsbench4.inp");
      std::map<int,BoundaryFunction<2>* > shapes;
      shapes[80] = &RK;

      BasicInit(inpname,prerefine,shapes);
    }
  ~BenchMarkMeshAgent() {}
};

/*---------------------------------------------------*/

#endif
