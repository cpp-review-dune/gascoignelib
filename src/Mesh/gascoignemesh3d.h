#ifndef  __GascoigneMesh3d_h
#define  __GascoigneMesh3d_h

#include  "gascoignemesh.h"
#include  "vertex.h"

/*-----------------------------------------*/

class GascoigneMesh3d : public GascoigneMesh
{
protected:

  // basic
  std::vector<Vertex3d>   nx;

public:

  GascoigneMesh3d();
  ~GascoigneMesh3d() {}

  std::vector<Vertex3d>& GetVertexVector() {return nx;}

  int  dimension() const {return 3;}
  int  nnodes()    const {return nx.size();}
  int  ncells()    const {return nc.size()/8;}
  int  nodes_per_cell()  const {return 8;}
  const Vertex3d& vertex3d(int i) const { return nx[i];} 
  int  vertex_of_cell(int i, int ii) const { return nc[8*i+ii]; }

  nvector<int>  IndicesOfCell(int iq) const;
};

#endif
