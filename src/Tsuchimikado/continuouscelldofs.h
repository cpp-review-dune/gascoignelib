/*----------------------------   continuouscelldofs.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __continuouscelldofs_H
#define __continuouscelldofs_H
/*----------------------------   continuouscelldofs.h     ---------------------------*/

#include "celldofs.h"
#include "triacontainer.h"
#include "meshlevel.h"

namespace Tsuchimikado
{

  /**
   * continuous degrees of freedom in cells
   * Template: DIM dimension of the domain 
   * and PD patchdepth, number of global refinements
   **/

  template<int DIM, int PD>
    class ContinuousCellDofs : public CellDofs<DIM>
  {
    const TriaContainer<DIM>* __TC;
    
  public:
    ContinuousCellDofs()
      {
	__TC = 0;
      }
        
    void BasicInit(const TriaContainer<DIM>* TC)
    {
      __TC = TC;
    }
    
    void ReInit(const MeshLevel<DIM>& ML);
    
  };
  


}



/*----------------------------   continuouscelldofs.h     ---------------------------*/
/* end of #ifndef __continuouscelldofs_H */
#endif
/*----------------------------   continuouscelldofs.h     ---------------------------*/
