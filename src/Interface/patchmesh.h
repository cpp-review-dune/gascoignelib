#ifndef  __PatchMesh_h
#define  __PatchMesh_h

#include  "meshinterface.h"

/*-----------------------------------------*/

class PatchMesh : public MeshInterface
{
  
 public:
  
  PatchMesh() {};
  virtual ~PatchMesh() {}
  
  virtual bool       HasPatch()                        const { assert(0);}
  virtual int        npatches()                          const { return 0; }
  virtual int        nodes_per_patch()                 const { assert(0);}
  virtual Gascoigne::IntVector  CoarseIndices(int iq)             const { assert(0);}
  virtual const Gascoigne::IntVector& IndicesOfPatch(int i)         const { assert(0);}
  virtual const Gascoigne::IntVector& VertexOnBoundary(int color) const { assert(0);}
  virtual const Gascoigne::IntVector& CellOnBoundary(int color)   const { assert(0);}
  virtual const Gascoigne::IntVector& LocalOnBoundary(int color)  const { assert(0);}

  virtual bool CellIsCurved(int iq)                    const { return 0;}

    // MPI
  virtual void send(int p) const { assert(0); }  
  virtual void recv(int p)       { assert(0); }
};

#endif
