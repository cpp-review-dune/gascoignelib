#ifndef  __GascoigneMesh2d_h
#define  __GascoigneMesh2d_h

#include  "gascoignemesh.h"
#include  "vertex.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMesh2d : public GascoigneMesh
{
protected:
  
  // basic
  std::vector<Vertex2d>   nx;
  
public:
  
  GascoigneMesh2d();
  ~GascoigneMesh2d();

  std::string GetName() const {return "GascoigneMesh2d";}

  std::vector<Vertex2d>& GetVertexVector() {return nx;}

  int  dimension() const {return 2;}
  int  nnodes()    const {return nx.size();}
  int  ncells()    const {return nc.size()/4;}
  int  nhanging()  const { return HangingHandler.GetStructure()->size(); }

  int  nodes_per_cell(int i)  const { return 4;}
  int  VtkType(int i) const { return 9;}

  const Vertex2d& vertex2d(int i) const { return nx[i];}
  int  vertex_of_cell(int i, int ii) const { return nc[4*i+ii]; }
  
  IntVector  IndicesOfCell(int iq) const;
};
}

#endif
