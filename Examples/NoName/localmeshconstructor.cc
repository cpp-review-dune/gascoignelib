#include  "localmeshconstructor.h"
#include  "localmesh2d.h"

/*-----------------------------------------*/

LocalMeshConstructor::LocalMeshConstructor(const HierarchicalMesh* mm,
					   GascoigneMultiGridMesh* gmg)
  : GascoigneMeshConstructor(mm,gmg)
{
}

/*-----------------------------------------*/

void LocalMeshConstructor::Construct2d(GascoigneMesh* NM, const LevelMesh2d* LM) const
{
  GascoigneMeshConstructor::Construct2d(NM,LM);
  LocalMesh2d* MP = dynamic_cast<LocalMesh2d*>(NM);
  assert(MP);
  MP->SetCoordinates(LM);
}

