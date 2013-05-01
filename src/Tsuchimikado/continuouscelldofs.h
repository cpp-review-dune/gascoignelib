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


  // dofs per line, quad, hex and element
#define DPL(p) ( (1<<p)+1 )
#define DPQ(p) ( ((1<<p)+1)*((1<<p)+1) )
#define DPH(p) ( ((1<<p)+1)*((1<<p)+1)*((1<<p)+1) )
#define DPE(d,p) ( (d==2)?(DPQ(p)):(DPH(p)) )
  
  template<int DIM, int PD>
    class ContinuousCellDofs : public CellDofs<DIM>
  {
    const TriaContainer<DIM>* __TC;

    // 
    std::vector<Gascoigne::fixarray<DPE(DIM,PD), int> > __dofs; 

    typedef std::tr1::unordered_map<int,int> HASH_MAP;
    typedef std::tr1::unordered_set<int>     HASH_SET;
    
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
    void create_dofs_in_elements_2d(const MeshLevel<DIM>& ML, const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
				   const std::vector<int>& idol, const std::vector<int>& idoq);
    void create_dofs_in_elements_3d(const MeshLevel<DIM>& ML, const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
				    const std::vector<int>& idol, const std::vector<int>& idoq, const std::vector<int>& idoh);
    std::vector<int> create_dofs_on_one_line(int cell, int li, const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
					     const std::vector<int>& idol);
    

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
