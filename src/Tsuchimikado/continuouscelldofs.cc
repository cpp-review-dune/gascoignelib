#include "continuouscelldofs.h"
#include <iostream>

using namespace std;
using namespace Gascoigne;


namespace Tsuchimikado
{

  template<int DIM, int PD>
  void ContinuousCellDofs<DIM,PD>::ReInit(const MeshLevel<DIM>& ML)
  {
    assert((DIM==2)||(DIM==3));
    this->__vertices.clear();
    
    // number of vertices in a line
    int nvl = (1<<PD)+1;
    // number of vertices in a cell
    int nvc = nvl*nvl; if (DIM==3) nvc*=nvl;
    
    
    // create dofs & vertices
    for (int c=0;c<ML.size();++c)
      {
	if (DIM==2)      create_dofs_in_element2d(ML[c]);
	else if (DIM==3) create_dofs_in_element3d(ML[c]);
	else abort();
      }
  }


  template<int DIM, int PD>
  Vertex<DIM> ContinuousCellDofs<DIM,PD>::create_dofs_in_element2d(int cell)
  {

  }
  

  template<int DIM, int PD>
  Vertex<DIM> ContinuousCellDofs<DIM,PD>::create_dofs_in_element3d(int cell)
  {
  }
  





  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  template<int DIM, int PD>
  Vertex<DIM> ContinuousCellDofs<DIM,PD>::CreateVertex2d(int cell, int ix,int iy) const
  {
    assert(DIM==2);
    assert(__TC);
    assert(cell<__TC->ncells());
    
    Vertex<DIM> v;
    double sx = static_cast<double>(ix)/static_cast<double> (dofs_per_line()-1);
    double sy = static_cast<double>(iy)/static_cast<double> (dofs_per_line()-1);
    assert((sx>=0.0)&&(sx<=1.0));
    assert((sy>=0.0)&&(sy<=1.0));
    v.equ((1.0-sx)*(1.0-sy) , __TC->vertex(__TC->cell(cell).node(0)),
	  sx      *(1.0-sy) , __TC->vertex(__TC->cell(cell).node(1)), 
	  sx      *sy       , __TC->vertex(__TC->cell(cell).node(2)), 
	  (1.0-sx)*sy       , __TC->vertex(__TC->cell(cell).node(3)));
    return v;
  }

  template<int DIM, int PD>
  Vertex<DIM> ContinuousCellDofs<DIM,PD>::CreateVertex3d(int cell, int ix,int iy, int iz) const

  {
    assert(DIM==3);
    assert(__TC);
    assert(cell<__TC->ncells());
    
    Vertex<DIM> v;
    double sx = static_cast<double>(ix)/static_cast<double> (dofs_per_line()-1);
    double sy = static_cast<double>(iy)/static_cast<double> (dofs_per_line()-1);
    double sz = static_cast<double>(iz)/static_cast<double> (dofs_per_line()-1);
    assert((sx>=0.0)&&(sx<=1.0));
    assert((sy>=0.0)&&(sy<=1.0));
    assert((sz>=0.0)&&(sz<=1.0));
    v.equ((1.0-sx)*(1.0-sy)*(1.0-sz) , __TC->vertex(__TC->cell(cell).node(0)),
	  sx      *(1.0-sy)*(1.0-sz) , __TC->vertex(__TC->cell(cell).node(1)), 
	  sx      *sy      *(1.0-sz) , __TC->vertex(__TC->cell(cell).node(2)), 
	  (1.0-sx)*sy      *(1.0-sz) , __TC->vertex(__TC->cell(cell).node(3)),
	  (1.0-sx)*(1.0-sy)*sz       , __TC->vertex(__TC->cell(cell).node(4)),
	  sx      *(1.0-sy)*sz       , __TC->vertex(__TC->cell(cell).node(5)), 
	  sx      *sy      *sz       , __TC->vertex(__TC->cell(cell).node(6)), 
	  (1.0-sx)*sy      *sz       , __TC->vertex(__TC->cell(cell).node(7)));
    return v;
  }
  



  template class ContinuousCellDofs<2,0>;
  template class ContinuousCellDofs<2,1>;
  template class ContinuousCellDofs<2,2>;
  


}


