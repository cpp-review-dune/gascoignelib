#ifndef  __LocalMultiGridMesh_h
#define  __LocalMultiGridMesh_h

#include  "gascoignemultigridmesh.h"
#include  "localmesh2d.h"

/*-----------------------------------------*/

class LocalMultiGridMesh : public GascoigneMultiGridMesh
{
private:

  virtual GascoigneMesh* NewMesh(int dim) 
    {
      if     (dim==2) return new LocalMesh2d;
      assert(0);
    }

public:

  LocalMultiGridMesh() {}

};

#endif






