/*----------------------------   dgdofs.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgdofs_H
#define __dgdofs_H
/*----------------------------   dgdofs.h     ---------------------------*/

#include "meshinterface.h"
#include "nvector.h"
#include "dgbase.h"
#include "matrixinterface.h"

namespace Gascoigne
{
  
  template <class BASE>
  class DGDofHandler
  {
    typedef std::array<size_t, BASE::N> ElementType;
    /**
     * 
     * For every edge we store the two adjacent elements e1,e2
     * and the local orientation index within the element l1,l2
     * So, (e1,e2,l1,l2) for an internal edge and (e1,-1,l1,-1) 
     * for a boundary edge
     * The normal vector will always be computed based on e1
     *
     *  2D: (local index given on the edges)
     * 
     *  2 - 2 - 3
     *  |       |
     *  3       1
     *  |       |
     *  0 - 0 - 1
     *
     **/
    typedef std::array<size_t, 4>       EdgeType;
    const static size_t NONEDGE = -1;

  protected:
    
    ///// CellDofs
    size_t _ndofs;

    std::vector<ElementType> _elements;
    std::vector<EdgeType>    _edges;

  public:
    // INIT
    void InitFromGascoigneMesh(const MeshInterface *M);

    // ACCESS
    int nelements() const
    {
      return _elements.size();
    }
    int nedges() const
    {
      return _edges.size();
    }
    int ndofs() const
    {
      return _ndofs;
    }
    const ElementType &getelement(int e) const
    {
      assert(e < _elements.size());
      return _elements[e];
    }
    const EdgeType &getedge(int e) const
    {
      assert(e < _edges.size());
      return _edges[e];
    }

    //// Matrix Structure
    void Structure(SparseStructureInterface *S) const;

    ///// Global/Local
    void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const;
    void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const;
    void LocalToGlobalMatrix(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
    void LocalToGlobalMatrix(MatrixInterface& A, EntryMatrix& E, int iq1, int iq2, double s) const;
  };
}


/*----------------------------   dgdofs.h     ---------------------------*/
/* end of #ifndef __dgdofs_H */
#endif
/*----------------------------   dgdofs.h     ---------------------------*/
