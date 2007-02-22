#ifndef  __MeshInterface_h
#define  __MeshInterface_h

#include  "vertex.h"
#include  <set>
#include  <string>
#include  "gascoigne.h"
#include  "paramfile.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments MeshInterface

  ///
  ///
  /////////////////////////////////////////////

  class MeshInterface
  {
    private:

    protected:

    public:
      MeshInterface() {}
      virtual ~MeshInterface() {}

      virtual void BasicInit(const ParamFile* pf)=0;

      virtual std::string GetName() const=0;

      virtual int  dimension() const=0;
      virtual int  nnodes()    const=0;
      virtual int  nhanging()  const {
        return 0;
      }
      virtual int  ncells()    const=0;

      virtual int  nodes_per_cell(int i)         const=0;
      virtual int  vertex_of_cell(int i, int ii) const=0;
      virtual const Vertex2d& vertex2d(int i)    const {
        std::cerr << "\"MeshInterface::vertex2d\" not written!" << std::endl;
        abort();
      }
      virtual const Vertex3d& vertex3d(int i)    const {
        std::cerr << "\"MeshInterface::vertex3d\" not written!" << std::endl;
        abort();
      }
      virtual       IntSet     GetColors()                const=0;
      virtual const IntVector* Vertexo2n()                const=0;
      virtual       IntVector  IndicesOfCell(int iq)      const {
        std::cerr << "\"MeshInterface::IndicesOfCell\" not written!" << std::endl;
        abort();
      }
      virtual const IntVector* CellOnBoundary(int color)  const {
        std::cerr << "\"MeshInterface::CellOnBoundary\" not written!" << std::endl;
        abort();
      }
      virtual const IntVector* LocalOnBoundary(int color) const {
        std::cerr << "\"MeshInterface::LocalOnBoundary\" not written!" << std::endl;
        abort();
      }
      virtual const IntVector* PatchOnBoundary(int color)  const {
        std::cerr << "\"MeshInterface::PatchOnBoundary\" not written!" << std::endl;
        abort();
      }
      virtual const IntVector* LocalPatchOnBoundary(int color) const {
        std::cerr << "\"MeshInterface::LocalPatchOnBoundary\" not written!" << std::endl;
        abort();
      }
      virtual const IntVector* VertexOnBoundary(int col) const {
        std::cerr << "\"MeshInterface::VertexOnBoundary\" not written!" << std::endl;
        abort();
      }
      virtual int VtkType(int i) const=0;
  };
}

#endif
