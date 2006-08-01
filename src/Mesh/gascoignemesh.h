#ifndef  __GascoigneMesh_h
#define  __GascoigneMesh_h

#include  "patchmesh.h"
#include  "patchindexhandler.h"
#include  "boundaryindexhandler.h"
#include  "hangingindexhandler.h"
#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMesh : public PatchMesh
{
protected:

  IntVector              nc, vo2n;
  PatchIndexHandler      PatchHandler;
  BoundaryIndexHandler   BoundaryHandler;
  HangingIndexHandler    HangingHandler;

public:

  GascoigneMesh();
  virtual ~GascoigneMesh();

	void BasicInit(const ParamFile* pf) {assert(0);}

  const IntVector&               GetCellVector()          const  { return nc;}
  const PatchIndexHandler&    GetPatchIndexHandler()   const  { return PatchHandler;}
  const BoundaryIndexHandler&    GetBoundaryIndexHandler()const  { return BoundaryHandler;}
  const HangingIndexHandler&  GetHangingIndexHandler() const  { return HangingHandler;}
  const IntVector*               Vertexo2n()              const  { return &vo2n;}

  IntVector&               GetCellVector()           { return nc;}
  PatchIndexHandler&    GetPatchIndexHandler()    { return PatchHandler;}
  BoundaryIndexHandler&    GetBoundaryIndexHandler() { return BoundaryHandler;}
  HangingIndexHandler&  GetHangingIndexHandler()  { return HangingHandler;}
  IntVector*               Vertexo2n()               { return &vo2n;}

  int  nodes_per_patch()   const { return PatchHandler.nodes_per_patch();}
  int  nodes_per_q4patch() const { return PatchHandler.nodes_per_q4patch();}
  bool HasPatch()          const { return PatchHandler.HasPatch();}
  bool HasQ4Patch()        const { return PatchHandler.HasQ4Patch();}
  int  npatches()          const { return PatchHandler.npatches(); }
  int  nq4patches()        const { return PatchHandler.nq4patches(); }

  const IntVector* IndicesOfPatch    (int i)     const { return &(PatchHandler.IndicesOfPatch(i));}
  const IntVector* IndicesOfQ4Patch  (int i)     const { return &(PatchHandler.IndicesOfQ4Patch(i));}
  const IntVector* VertexOnBoundary(int color) const { return &(BoundaryHandler.Verteces(color)); }
  const IntVector* CellOnBoundary  (int color) const { return &(BoundaryHandler.Cells(color));    }
  const IntVector* LocalOnBoundary (int color) const { return &(BoundaryHandler.Localind(color)); }
  
  std::set<int> GetColors()             const { return BoundaryHandler.GetColors();}
  IntVector  Q2IndicesOfQ4Patch(int i)  const { return PatchHandler.Q2IndicesOfQ4Patch(i);}
  IntVector  CoarseIndices(int iq)   const { return PatchHandler.CoarseIndices(iq);}
  IntVector  CoarseIndicesQ4(int iq) const { return PatchHandler.CoarseIndicesQ4(iq);}

  // wird von GascoigneMesh2d/3d ueberschrieben
  virtual IntVector  IndicesOfCell(int iq) const{ assert(0); return IntVector(); }
  
  virtual bool CellIsCurved(int iq) const { return 0;}
  virtual int  nedges()    const { return 0; }
};
}

#endif
