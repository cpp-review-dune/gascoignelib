#ifndef  __PatchMesh_h
#define  __PatchMesh_h

#include  "meshinterface.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class PatchMesh : public virtual MeshInterface
  {
  
    private:

    protected:

    public:
      PatchMesh() {};
      virtual ~PatchMesh() {}
 
      virtual bool       HasPatch()                        const=0;
      virtual int        npatches()                        const=0;
      virtual int        nodes_per_patch()                 const=0;
      virtual IntVector  CoarseIndices(int iq)             const=0;
      virtual const IntVector* IndicesOfPatch(int i)       const=0;
      virtual const IntVector* VertexOnBoundary(int color) const=0;
      virtual const IntVector* CellOnBoundary(int color)   const=0;
      virtual const IntVector* LocalOnBoundary(int color)  const=0;

      virtual bool CellIsCurved(int iq)                    const {
        return false;
      }

      // MPI
      virtual void send(int p) const {
        std::cerr << "\"PatchMesh::send\" not written!" << std::endl;
        abort();
      }
      virtual void recv(int p) {
        std::cerr << "\"PatchMesh::recv\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
