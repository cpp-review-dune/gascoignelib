#ifndef  __LocalMeshConstructor_h
#define  __LocalMeshConstructor_h

#include  "gascoignemeshconstructor.h"
#include  "levelmesh2d.h"

/*-----------------------------------------*/

class LocalMeshConstructor : public GascoigneMeshConstructor
{
protected:

  void Construct2d(GascoigneMesh* NM, const LevelMesh2d* LM) const;

public:

  LocalMeshConstructor(const HierarchicalMesh* mm, GascoigneMultiGridMesh* gmg);
  
};

#endif






