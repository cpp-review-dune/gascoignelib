#ifndef __LocalMeshAgent_h
#define __LocalMeshAgent_h

#include "runderkreis.h"
#include "meshagent.h"

/*---------------------------------------------------*/

class LocalMeshAgent : public MeshAgent
{
 protected:
  
  RunderKreis RK;

 public:
  
  LocalMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      int dim=2;
      int prerefine=3;
      std::string inpname("../NavierStokes2D/nsbench4.inp");

      AddShape(80,&RK);
      BasicInit(dim,inpname,prerefine);
    }
  ~LocalMeshAgent() {}
};

/*---------------------------------------------------*/

#endif
