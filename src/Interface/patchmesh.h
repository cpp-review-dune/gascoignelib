#ifndef  __PatchMesh_h
#define  __PatchMesh_h

#include  "meshinterface.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PatchMesh : public MeshInterface
{
  
 public:
  
  PatchMesh() {};
  virtual ~PatchMesh() {}
  
  virtual bool       HasPatch()                        const { assert(0); return 0;}
  virtual int        npatches()                          const { return 0;}
  virtual int        nodes_per_patch()                 const { assert(0); return -1;}
  virtual IntVector  CoarseIndices(int iq)             const { assert(0); return IntVector();}
  virtual const IntVector* IndicesOfPatch(int i)         const { assert(0); return NULL;}
  virtual const IntVector* VertexOnBoundary(int color) const { assert(0); return NULL;}
  virtual const IntVector* CellOnBoundary(int color)   const { assert(0); return NULL;}
  virtual const IntVector* LocalOnBoundary(int color)  const { assert(0); return NULL;}

  virtual bool CellIsCurved(int iq)                    const { return 0;}

    // MPI
  virtual void send(int p) const { assert(0); }  
  virtual void recv(int p)       { assert(0); }
};
}

#endif
