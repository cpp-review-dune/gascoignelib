#include "continuouscelldofs.h"
#include <iostream>

using namespace std;
using namespace Gascoigne;


namespace Tsuchimikado
{

  template<int DIM, int PS>
  void ContinuousCellDofs<DIM,PS>::ReInit(const MeshLevel<DIM>& ML)
  {
    assert((DIM==2)||(DIM==3));
    this->__vertices.clear();
    
    // number of vertices in a line
    int nvl = (1<<PS)+1;
    // number of vertices in a cell
    int nvc = nvl*nvl; if (DIM==3) nvc*=nvl;
    
     
    
  }



  template class ContinuousCellDofs<2,0>;
  template class ContinuousCellDofs<2,1>;
  template class ContinuousCellDofs<2,2>;
  


}


