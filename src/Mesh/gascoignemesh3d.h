#ifndef  __GascoigneMesh3d_h
#define  __GascoigneMesh3d_h

#include  "gascoignemesh.h"
#include  "vertex.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMesh3d : public GascoigneMesh
{
protected:

  // basic
  std::vector<Vertex3d>   nx;

public:

  GascoigneMesh3d();
  ~GascoigneMesh3d() {}

  std::string GetName() const {return "GascoigneMesh3d";}

        std::vector<Vertex3d>& GetVertexVector()       {return nx;}
  const std::vector<Vertex3d>& GetVertexVector() const {return nx;}

  int  dimension() const {return 3;}
  int  nnodes()    const {return nx.size();}
  int  ncells()    const {return nc.size()/8;}

  int  nodes_per_cell(int i)  const { return 8;}
  int  VtkType(int i) const { return 12;}

  const Vertex3d& vertex3d(int i) const { return nx[i];} 
  int  vertex_of_cell(int i, int ii) const { return nc[8*i+ii]; }

  IntVector  IndicesOfCell(int iq) const;
};
}

#endif
