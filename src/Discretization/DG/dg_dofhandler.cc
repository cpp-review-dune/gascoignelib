#include "dg_dofhandler.h"

namespace Gascoigne
{
  
  // // get the dof-indices for one edge
  // template<int DIM>
  // void DGDofHandler<DIM>::GetIndicesEdge(int e, std::vector<int>& indices) const
  // {
  //   assert(indices.size()==this->ndofs_per_edge());
  //   assert(e < nedges());
    
  //   const DGEdge& E = this->edge(e);
  //   // get the two elements
  //   int c1 = E.master(); assert(c1!=-1);
  //   int c2 = E.slave();

  //   int l1 = E.localmaster();
  //   int l2 = E.localslave();

  //   abort();    
  // }

  template class DGDofHandler<2>;
  template class DGDofHandler<3>;

}
