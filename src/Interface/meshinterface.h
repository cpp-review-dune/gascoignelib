#ifndef  __MeshInterface_h
#define  __MeshInterface_h


/////////////////////////////////////////////
///
///@brief
///  ... comments MeshInterface

///
///
/////////////////////////////////////////////

#include  "vertex.h"
#include  <set>
#include  "gascoigne.h"

using namespace Gascoigne;

class MeshInterface
{
public:

  MeshInterface() {}
  virtual ~MeshInterface() {}

  virtual int  dimension() const=0;
  virtual int  nnodes()    const=0;
  virtual int  ncells()    const=0;

  virtual int  nodes_per_cell()         const { assert(0);}
  virtual int  vertex_of_cell(int, int)    const=0;
  virtual const Vertex2d& vertex2d(int i)  const { assert(0);}
  virtual const Vertex3d& vertex3d(int i)  const { assert(0);} 
  virtual std::set<int> GetColors()        const { assert(0);}
  virtual IntVector  IndicesOfCell(int iq) const { assert(0);}
  virtual const IntVector& Vertexo2n()     const { assert(0);}

  virtual const IntVector& CellOnBoundary(int color)   const { assert(0);}
  virtual const IntVector& LocalOnBoundary(int color)  const { assert(0);}
};


#endif
