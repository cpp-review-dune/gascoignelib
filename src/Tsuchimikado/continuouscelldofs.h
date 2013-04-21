/*----------------------------   continuouscelldofs.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __continuouscelldofs_H
#define __continuouscelldofs_H
/*----------------------------   continuouscelldofs.h     ---------------------------*/

#include "celldofs.h"
#include "triacontainer.h"
#include "meshlevel.h"
#include "vertex.h"

namespace Tsuchimikado
{

  /**
   * continuous degrees of freedom in cells
   * Template: DIM dimension of the domain 
   * and PD patchdepth, number of global refinements
   *
   * PD=0: Q1
   * PD=1: Q2
   * PD=2: Q4
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


    /// simple access
    int dofs_per_line() const
    { return (1<<PD)+1; }

    int dofs_per_element() const
    {
      if (DIM==2)      return dofs_per_line()*dofs_per_line();
      else if (DIM==3) return dofs_per_line()*dofs_per_line()*dofs_per_line();
      abort();
    }

    

  protected:
    //
    void create_dofs_in_element_2d(int cell);
    void create_dofs_in_element_3d(int cell);

    Gascoigne::Vertex<DIM> CreateVertex2d(int cell, int ix,int iy) const;
    Gascoigne::Vertex<DIM> CreateVertex3d(int cell, int ix,int iy, int iz) const;
    
  };
  


}



/*----------------------------   continuouscelldofs.h     ---------------------------*/
/* end of #ifndef __continuouscelldofs_H */
#endif
/*----------------------------   continuouscelldofs.h     ---------------------------*/
