#include "continuouscelldofs.h"
#include <iostream>

using namespace std;
using namespace Gascoigne;

#define DPL(p) ( (1<<p)+1 )
#define DPQ(p) ( ((1<<p)+1)*((1<<p)+1) )
#define DPH(p) ( ((1<<p)+1)*((1<<p)+1)*((1<<p)+1) )
#define DPE(d,p) ( (d==2)?(DPQ(p)):(DPH(p)) )


namespace Tsuchimikado
{

  template<int DIM, int PD>
  void ContinuousCellDofs<DIM,PD>::ReInit(const MeshLevel<DIM>& ML)
  {
    assert((DIM==2)||(DIM==3));
    this->__vertices.clear();
    this->__dofs.clear();

    /////////// create dofs on nodes
    // list of active nodes
    int newdofs = ML.nnodes();
    // node g2l
    HASH_MAP  nodeg2l;
    int nn=0;
    for (typename MeshLevel<DIM>::CIT cit = ML.nodes_begin();cit!=ML.nodes_end();++cit,++nn)
      nodeg2l[*cit]=nn;
    assert(nn=ML.nnodes());

    HASH_MAP  lineg2l;
    nn=0;
    for (typename MeshLevel<DIM>::CIT cit = ML.lines_begin();cit!=ML.lines_end();++cit,++nn)
      lineg2l[*cit]=nn;
    assert(nn=ML.nlines());

    /////////// create dofs on lines
    std::vector<int> inner_dofs_on_lines;
    if (PD>0)
      {
	for (std::vector<int>::const_iterator it = ML.lines_begin();
	     it!=ML.lines_end();++it)
	  {
	    if (ML.line_is_hanging(*it))
	      for (int i=0;i<DPL(PD)-2;++i)
		inner_dofs_on_lines.push_back(-1);
	    else
	      for (int i=0;i<DPL(PD)-2;++i)
		inner_dofs_on_lines.push_back(newdofs++);
	  }

	// one hanging lines, set inner dofs
	const HASH_SET& HL = ML.GetHangingLines();
	for (HASH_SET::const_iterator hlit = HL.begin(); hlit != HL.end();++hlit)
	  {
	    assert(lineg2l.find(*hlit)!=lineg2l.end());
	    int loc = lineg2l[*hlit];
	    assert(__TC->line(*hlit).nchilds()==2);

	    int lc0 = __TC->line(*hlit).child(0); assert(lineg2l.find(lc0)!=lineg2l.end());
	    int loc0 = lineg2l[lc0];
	    int lc1 = __TC->line(*hlit).child(1); assert(lineg2l.find(lc1)!=lineg2l.end());
	    int loc1 = lineg2l[lc1];
	    //
	    for (int i=1;i<(DPL(PD)-1)/2;++i) // left half
	      {
		assert(inner_dofs_on_lines[(DPL(PD)-2)*loc+i-1]==-1);
		assert(inner_dofs_on_lines[(DPL(PD)-2)*loc0+2*i-1]!=-1);
		inner_dofs_on_lines[(DPL(PD)-2)*loc+i-1]=inner_dofs_on_lines[(DPL(PD)-2)*loc0+2*i-1];
	      }
	    for (int i=1;i<(DPL(PD)-1)/2;++i) // right half
	      {
		assert(inner_dofs_on_lines[(DPL(PD)-2)*loc+i-1+ (DPL(PD)-1)/2 ]==-1);
		assert(inner_dofs_on_lines[(DPL(PD)-2)*loc1+2*i-1]!=-1);
		inner_dofs_on_lines[(DPL(PD)-2)*loc+i-1 + (DPL(PD)-1)/2]=inner_dofs_on_lines[(DPL(PD)-2)*loc1+2*i-1];
	      }
	    // middle
	    assert(inner_dofs_on_lines[(DPL(PD)-2)*loc + (DPL(PD)-1)/2-1]==-1);
	    int mn = __TC->line(lc0).node(1);
	    assert(nodeg2l.find(mn)!=nodeg2l.end());
	    inner_dofs_on_lines[(DPL(PD)-2)*loc + (DPL(PD)-1)/2-1] = nodeg2l[mn];
	  }
      }
    
    
    

    /////////// create dofs on quads
    std::vector<int> inner_dofs_on_quads;
    if (PD>0)
      for (std::vector<int>::const_iterator it = ML.quads_begin();
	   it!=ML.quads_end();++it)
	{
	  if (ML.quad_is_hanging(*it))
	    for (int i=0;i<(DPL(PD)-2)*(DPL(PD)-2);++i)
	      inner_dofs_on_quads.push_back(-1);
	  else
	    for (int i=0;i<(DPL(PD)-2)*(DPL(PD)-2);++i)
	      inner_dofs_on_quads.push_back(newdofs++);
	}

    /////////// create dofs on hexes
    std::vector<int> inner_dofs_on_hexes;
    if (PD>0)
      for (std::vector<int>::const_iterator it = ML.hexes_begin();
	   it!=ML.hexes_end();++it)
	for (int i=0;i<(DPL(PD)-2)*(DPL(PD)-2)*(DPL(PD)-2);++i)
	  inner_dofs_on_hexes.push_back(newdofs++);

    // create dofs & vertices
    if (DIM==2)
      create_dofs_in_elements_2d(ML,nodeg2l,lineg2l,inner_dofs_on_lines, inner_dofs_on_quads);
    else if (DIM==3)
      create_dofs_in_elements_3d(ML,nodeg2l,lineg2l,inner_dofs_on_lines, inner_dofs_on_quads, inner_dofs_on_hexes);
    else abort();


    for (int i=0;i<__dofs.size();++i)
      for (int j=0;j<__dofs[i].size();++j)
	assert(__dofs[i][j]>=0);
    
  }
  

  template<int DIM, int PD>
  void ContinuousCellDofs<DIM,PD>::create_dofs_in_elements_2d(const MeshLevel<DIM>& ML, const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
							      const vector<int>& idol, const vector<int>& idoq)
  {
    assert(DIM==2);
    __dofs.clear();
    // create lexicographic ordering of dofs
    int nquad=0;
    for (typename MeshLevel<DIM>::CIT cit = ML.begin();cit!=ML.end();++cit,++nquad)
      {
	fixarray<DPE(DIM,PD), int> dof;

	///////////////// dofs on outer nodes
	HASH_MAP::const_iterator nit; 
	nit = nodeg2l.find(__TC->cell(*cit).node(0)); assert(nit!=nodeg2l.end());
	dof[0]                   = nit->second;

	nit = nodeg2l.find(__TC->cell(*cit).node(1)); assert(nit!=nodeg2l.end());
	dof[DPL(PD)-1]           = nit->second;

	nit = nodeg2l.find(__TC->cell(*cit).node(3)); assert(nit!=nodeg2l.end());
	dof[DPL(PD)*(DPL(PD)-1)] = nit->second;

	nit = nodeg2l.find(__TC->cell(*cit).node(2)); assert(nit!=nodeg2l.end());
	dof[DPL(PD)*DPL(PD)-1]   = nit->second;

	
	////////////////// dofs on the lines
	vector<int> dool;
	dool=create_dofs_on_one_line(*cit,0,nodeg2l,lineg2l,idol);
	for (int ix=1;ix<DPL(PD)-1;++ix) 
	  dof[ix] = dool[ix-1];

	dool=create_dofs_on_one_line(*cit,1,nodeg2l,lineg2l,idol);
	for (int iy=1;iy<DPL(PD)-1;++iy) 
	  dof[iy*DPL(PD)+DPL(PD)-1] = dool[iy-1];

	dool=create_dofs_on_one_line(*cit,2,nodeg2l,lineg2l,idol);
	for (int ix=1;ix<DPL(PD)-1;++ix) 
	  dof[DPL(PD)*DPL(PD)-1-ix] = dool[ix-1];

	dool=create_dofs_on_one_line(*cit,3,nodeg2l,lineg2l,idol);
	for (int iy=1;iy<DPL(PD)-1;++iy) 
	  dof[(DPL(PD)-1-iy)*DPL(PD)] = dool[iy-1];
	




	
	////////////////// inner dofs
	for (int iy=1;iy<DPL(PD)-1;++iy)
	  for (int ix=1;ix<DPL(PD)-1;++ix)
	    dof[iy*DPL(PD)+ix] = idoq[nquad*(DPL(PD)-2)*(DPL(PD)-2) + (iy-1)*(DPL(PD)-2) + (ix-1) ];

	__dofs.push_back(dof);
      }
    
  }

  template<int DIM, int PD>
  void ContinuousCellDofs<DIM,PD>::create_dofs_in_elements_3d(const MeshLevel<DIM>& ML, const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
							      const vector<int>& idol, const vector<int>& idoq, const vector<int>& idoh)
  {
    assert(DIM==3);
    abort();
  }
  
  



  template<int DIM,int PD>
  std::vector<int> ContinuousCellDofs<DIM,PD>::create_dofs_on_one_line(int cell, int li, 
								       const HASH_MAP& nodeg2l, const HASH_MAP& lineg2l, 
								       const vector<int>& idol)
  {
    if (PD==0) return vector<int>();

    int line = __TC->cell(cell).line(li);
    HASH_MAP::const_iterator lineit = lineg2l.find(line);
    assert(lineit!=lineg2l.end());
    int localline = lineit->second;
    
    int n0 = __TC->line(line).node(0);
    int n1 = __TC->line(line).node(1);
    
    int dir=-1;
    if (__TC->cell(cell).node(li) == n0) dir = 0;
    if (__TC->cell(cell).node(li) == n1) dir = 1;
    assert(dir>=0);
    
    vector<int> dool;
    dool.resize(DPL(PD)-2,-1);

    if (dir==0)
      for (int ix=1;ix<DPL(PD)-1;++ix)
	{
	  assert(localline*(DPL(PD)-2)+ix-1<idol.size());
	  dool[ix-1] = idol[localline*(DPL(PD)-2)+ix-1];
	}
    else 
      for (int ix=1;ix<DPL(PD)-1;++ix)
	{
	  assert(localline*(DPL(PD)-2)+ix-1<idol.size());
	  dool[DPL(PD)-2-ix] = idol[localline*(DPL(PD)-2)+ix-1];
	}
    

    return dool;
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
  // template class ContinuousCellDofs<2,2>;
  // template class ContinuousCellDofs<2,3>;
  // template class ContinuousCellDofs<2,4>;
  


}



#undef DPL
#undef DPQ
#undef DPH
#undef DPE
