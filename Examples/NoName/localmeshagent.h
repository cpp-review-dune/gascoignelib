#ifndef  __LocalMeshAgent_h
#define  __LocalMeshAgent_h

#include  "meshagent.h"
#include  "localmeshconstructor.h"
#include  "localmultigridmesh.h"

/*-----------------------------------------*/

class LocalMeshAgent : public MeshAgent
{
 protected:
  
  GascoigneMultiGridMesh* NewMultiGridMesh() { return new LocalMultiGridMesh;}  
  
  void ReInit()
    {
      GMG->ReInit(dimension,HMP->nlevels()-1);
      LocalMeshConstructor MGM(HMP,GMG);
      MGM.BasicInit();
    }

 public:
  
  LocalMeshAgent() : MeshAgent() {}

  void refine_nodes(nvector<int>& refnodes, nvector<int>& coarsenodes) {
    std::cerr << "§§§\n";
    assert(HMP);
    HMP->vertex_patch_refine(refnodes,coarsenodes);
    MeshAgent::ReInit();
  }
  void SetCoordinates()
  {
    LocalMeshConstructor MGM(HMP,GMG);
    MGM.BasicInit();  
  }
};

#endif






