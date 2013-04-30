/*----------------------------   continuouscelldofs.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __continuouscelldofs_H
#define __continuouscelldofs_H
/*----------------------------   continuouscelldofs.h     ---------------------------*/

#include "celldofs.h"
#include "triacontainer.h"
#include "meshlevel.h"
#include "vertex.h"
#include "fixarray.h"
#include <vector>

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


  // dofs per line and dofs per element
#define DPL(p) ( (1<<p)+1 )
#define DPE(d,p) ( (d==2)?(DPL(p)*DPL(p)):(DPL(p)*DPL(p)*DPL(p)))
  
  template<int DIM, int PD>
    class ContinuousCellDofs : public CellDofs<DIM>
  {
    const TriaContainer<DIM>* __TC;

    std::vector<Gascoigne::fixarray<DPE(DIM,PD), int> > __dofs;
    
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
    int dofs_per_line()    const { return DPL(PD); }
    int dofs_per_element() const { return DPE(DIM,PD); }
    

  protected:
    //
    void create_dofs_in_element_2d(int cell);
    void create_dofs_in_element_3d(int cell);

    Gascoigne::Vertex<DIM> CreateVertex2d(int cell, int ix,int iy) const;
    Gascoigne::Vertex<DIM> CreateVertex3d(int cell, int ix,int iy, int iz) const;
    
  };
  


}

#undef DPL
#undef DPE


/*----------------------------   continuouscelldofs.h     ---------------------------*/
/* end of #ifndef __continuouscelldofs_H */
#endif
/*----------------------------   continuouscelldofs.h     ---------------------------*/
