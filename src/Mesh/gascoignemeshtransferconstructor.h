#ifndef  __GascoigneMeshTransferConstructor_h
#define  __GascoigneMeshTransferConstructor_h

#include  "gascoignemeshtransfer.h"
#include  "levelmesh2d.h"
#include  "levelmesh3d.h"

/*-----------------------------------------*/

class GascoigneMeshTransferConstructor2d
{
public:

  GascoigneMeshTransferConstructor2d(const HierarchicalMesh2d* HM, GascoigneMeshTransfer* GMT, 
				     const LevelMesh2d* LMfine, const LevelMesh2d* LMcoarse);
};

/*-----------------------------------------*/

class GascoigneMeshTransferConstructor3d
{
public:

  GascoigneMeshTransferConstructor3d(const HierarchicalMesh3d* HM, GascoigneMeshTransfer* GMT, 
				     const LevelMesh3d* LMfine, const LevelMesh3d* LMcoarse);
};

/*-----------------------------------------*/

#endif
