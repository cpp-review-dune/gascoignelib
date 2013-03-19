/*----------------------------   meshhierarchy.h     ---------------------------*/
/*      $Id: meshhierarchy.h,v 1.3 2008/11/04 15:47:17 meidner Exp $                 */
#ifndef __meshhierarchy_H
#define __meshhierarchy_H
/*----------------------------   meshhierarchy.h     ---------------------------*/


#include "meshlevel.h"
#include <vector>

namespace Tsuchimikado
{

  /**
   * The MeshHierarchy is a list of meshes for the multigrid-method.
   *
   * public functions:
   *
   *  ReInit(): initializes mesh hierarchy by global coarsening from TC
   *
   **/
  
  template<int DIM>
    class MeshHierarchy : public std::vector<MeshLevel<DIM> >
    {
    protected:
      
      const TriaContainer<DIM>&            __TC;
      
    public:
      
    MeshHierarchy(const TriaContainer<DIM>& TC) : __TC(TC){}
      
      void ReInit();
      int  nlevels() const { return this->size(); }
      
      MeshLevel<DIM>&       GetMeshLevel(int l)
	{ assert(l<this->size()); return std::vector<MeshLevel<DIM> >::operator[](l);}
      const MeshLevel<DIM>& GetMeshLevel(int l) const
      { assert(l<this->size()); return std::vector<MeshLevel<DIM> >::operator[](l);}

      void print_gnuplot(const std::string& fname) const;
      
    };
}



/*----------------------------   meshhierarchy.h     ---------------------------*/
/* end of #ifndef __meshhierarchy_H */
#endif
/*----------------------------   meshhierarchy.h     ---------------------------*/
