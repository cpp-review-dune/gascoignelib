#ifndef  __GascoigneMesh_h
#define  __GascoigneMesh_h

#include  "patchmesh.h"
#include  "patchindexhandler.h"
#include  "boundaryindexhandler.h"
#include  "hangingindexhandler.h"
#include  "gascoigne.h"

/*-----------------------------------------*/

class GascoigneMesh : public PatchMesh
{
protected:

  int                    mylevel;
  Gascoigne::IntVector              nc, vo2n;
  PatchIndexHandler      PatchHandler;
  BoundaryIndexHandler   BoundaryHandler;
  HangingIndexHandler    HangingHandler;

public:

  GascoigneMesh();
  virtual ~GascoigneMesh();

  const Gascoigne::IntVector&               GetCellVector()          const  { return nc;}
  const PatchIndexHandler&    GetPatchIndexHandler()   const  { return PatchHandler;}
  const BoundaryIndexHandler&    GetBoundaryIndexHandler()const  { return BoundaryHandler;}
  const HangingIndexHandler&  GetHangingIndexHandler() const  { return HangingHandler;}
  const Gascoigne::IntVector&               Vertexo2n()              const  { return vo2n;}

  Gascoigne::IntVector&               GetCellVector()           { return nc;}
  PatchIndexHandler&    GetPatchIndexHandler()    { return PatchHandler;}
  BoundaryIndexHandler&    GetBoundaryIndexHandler() { return BoundaryHandler;}
  HangingIndexHandler&  GetHangingIndexHandler()  { return HangingHandler;}
  Gascoigne::IntVector&               Vertexo2n()               { return vo2n;}

  int  nodes_per_patch() const { return PatchHandler.nodes_per_patch();}
  bool HasPatch()        const { return PatchHandler.HasPatch();}
  int  npatches()        const { return PatchHandler.npatches(); }

  const Gascoigne::IntVector& IndicesOfPatch    (int i)     const { return PatchHandler.IndicesOfPatch(i);}
  const Gascoigne::IntVector& VertexOnBoundary(int color) const { return BoundaryHandler.Verteces(color); }
  const Gascoigne::IntVector& CellOnBoundary  (int color) const { return BoundaryHandler.Cells(color);    }
  const Gascoigne::IntVector& LocalOnBoundary (int color) const { return BoundaryHandler.Localind(color); }
  
  std::set<int> GetColors()             const { return BoundaryHandler.GetColors();}
  Gascoigne::IntVector  CoarseIndices(int iq) const { return PatchHandler.CoarseIndices(iq);}

  // wird von GascoigneMesh2d/3d ueberschrieben
  virtual nvector<int>  IndicesOfCell(int iq) const{ assert(0);}
  
  virtual bool CellIsCurved(int iq) const { return 0;}
  virtual int  nedges()    const { return 0; }
};


#endif
