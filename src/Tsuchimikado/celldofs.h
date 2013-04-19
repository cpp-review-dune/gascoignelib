/*----------------------------   celldofs.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __celldofs_H
#define __celldofs_H
/*----------------------------   celldofs.h     ---------------------------*/


#include "vertex.h"
#include <vector>

namespace Tsuchimikado
{

  /**
   * Base class for the description of cell-wise degrees of freedom.
   * Used for all standard discretizations like Q1, Q2, .. continuous
   * or discontinuous, patched discretizations, LPS, and so on.
   **/

  template<int DIM>
    class CellDofs
    {
    protected:
      std::vector<Gascoigne::Vertex<DIM> > __vertices;

    public:

      
      int nvertices() const { return __vertices.size(); }
      
      const Gascoigne::Vertex<DIM>& vertex(int i) const { assert(i<nvertices()); return __vertices[i]; }
      Gascoigne::Vertex<DIM>& vertex(int i)             { assert(i<nvertices()); return __vertices[i]; }
      
    };
      
  
}


/*----------------------------   celldofs.h     ---------------------------*/
/* end of #ifndef __celldofs_H */
#endif
/*----------------------------   celldofs.h     ---------------------------*/
