#ifndef  __GascoigneMeshConstructor_h
#define  __GascoigneMeshConstructor_h

#include  "levelmesh2d.h"
#include  "levelmesh3d.h"
#include  "gascoignemultigridmesh.h"
#include  "gascoignemesh2d.h"
#include  "gascoignemesh3d.h"
#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMeshConstructor
{
protected:
  
  const HierarchicalMesh* HM;
  GascoigneMultiGridMesh* GMG;

  bool finestlevel;

  virtual void PatchToCell2d(PatchIndexHandler& PIH, const LevelMesh2d* LM) const;
  virtual void PatchToCell3d(PatchIndexHandler& PIH, const LevelMesh3d* LM) const;

  virtual void Construct2d(GascoigneMesh* NM, const LevelMesh2d* LM) const;
  virtual void Construct3d(GascoigneMesh* NM, const LevelMesh3d* LM) const;

  virtual LevelMesh2d* GascoigneMeshConstructor::LevelUpdate2d(GascoigneMesh* GM,
						     const IntSet& newquads, 
						     const IntSet& oldquads) const;
  virtual LevelMesh3d* GascoigneMeshConstructor::LevelUpdate3d(GascoigneMesh* GM,
						     const IntSet& newquads, 
						     const IntSet& oldquads) const;
  virtual void Loop2d();
  virtual void Loop3d();

public:
  
  GascoigneMeshConstructor(const HierarchicalMesh* mm,
			   GascoigneMultiGridMesh* gmg);

  virtual void BasicInit();
};
}

#endif
