#include <map>
#include "triacontainer.h"
#include "tsuchimikado.h"
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

#ifdef NDEBUG
#error assert must be enabled
#endif


namespace Tsuchimikado
{
  template<int DIM>
  void TriaContainer<DIM>::AddCurvedShape(int col, const BoundaryFunction<DIM>* BF)
  {
    assert(__curved_boundary.find(col)==__curved_boundary.end());
    __curved_boundary[col] = BF;
    std::cout << "cb: " << __curved_boundary.size() << std::endl;
  }
  
  template<int DIM>
  const std::map<int,const BoundaryFunction<DIM>* >& TriaContainer<DIM>::GetCurvedShapes() const
  {
    return __curved_boundary;
  }
  
  template<int DIM>
  bool TriaContainer<DIM>::isotropic() const
  { return (__mesh_flags&MESH_IS_ISOTROPIC); }

  template<int DIM>
  TriaContainer<DIM>::TriaContainer(unsigned int flags) : __mesh_flags(flags)
  {
  }

  //////////////////////////////
  ////////////////////////////// Init
  //////////////////////////////

  template<int DIM>
  void TriaContainer<DIM>::reset()
  {
    __vertices.clear();
    __hexes.clear();
    __quads.clear();
    __lines.clear();
    __hanging_lines.clear();
    __hanging_quads.clear();
    __boundary_lines.clear();
    __boundary_quads.clear();
  }
  

  //////////////////////////////
  ////////////////////////////// Access
  //////////////////////////////
  
  template<> const int TriaContainer<1>::ncells() const { abort(); }
  template<> const int TriaContainer<2>::ncells() const { return nquads(); }
  template<> const int TriaContainer<3>::ncells() const { return nhexes(); }


  template<> const TriaContainer<1>::CELL& TriaContainer<1>::cell(int i) const { abort(); }
  template<> const TriaContainer<2>::CELL& TriaContainer<2>::cell(int i) const { return quad(i); }
  template<> const TriaContainer<3>::CELL& TriaContainer<3>::cell(int i) const { return hex(i); }
  
  template<> TriaContainer<1>::CELL& TriaContainer<1>::cell(int i) { abort(); }
  template<> TriaContainer<2>::CELL& TriaContainer<2>::cell(int i) { return quad(i); }
  template<> TriaContainer<3>::CELL& TriaContainer<3>::cell(int i) { return hex(i); }

  // **************************************************
  
  template<int DIM>
  bool TriaContainer<DIM>::is_hanging(const LINE& L) const
  {
    if (L.nchilds()==0) return false;
    return node_is_hanging(middle_node(L));
  }
  template<int DIM>
  bool TriaContainer<DIM>::is_hanging(const QUAD& Q) const
  {
    if (Q.nchilds()==0) return false;
    return node_is_hanging(middle_node(Q));
  }
 
  // **************************************************

  
  // find an element
  template<int DIM>
  const int TriaContainer<DIM>::find(const HEX& H) const
  {
    for (int i=0;i<nhexes();++i) if (hex(i)==H) return i;
    return -1;
  }	  // **************************************************
  template<int DIM>
  const int TriaContainer<DIM>::find(const QUAD& Q) const
  {
    for (int i=0;i<nquads();++i) if (quad(i)==Q) return i;
    return -1;
  }	  // **************************************************
  template<int DIM>
  const int TriaContainer<DIM>::find(const LINE& L) const
  {
    for (int i=0;i<nlines();++i) if (line(i)==L) return i;
    return -1;
  }	  // **************************************************
  
  // middle nodes on lines, quads, hexes
  template<int DIM>
  const int TriaContainer<DIM>::middle_node(const HEX&  H) const
  {
    assert(H.type()==7);
    const HEX& C = hex(H.child(0));
    return C.node(6);
  }
  template<int DIM>
  const int TriaContainer<DIM>::middle_node(const QUAD& Q) const
  {
    // depends on type
    // for anisotropic quads this function crashes, if no middle-node exists
    assert(Q.type()>0);
   
    const QUAD& C = quad(Q.child(0));
    if      (Q.type()==3) return C.node(2);
    else if (Q.type()==1) return middle_node(line(C.line(1)));
    else if (Q.type()==2) return middle_node(line(C.line(2)));
    assert(0);
  }
  template<int DIM>
  const int TriaContainer<DIM>::middle_node(const LINE& L) const
  {
    assert(L.nchilds()==2);
    return line(L.child(0)).node(1);
  }

  // **************************************************

  
  template<int DIM>
  const int TriaContainer<DIM>::find_local_line_index_of_quad(int quad_number, int line_number) const
  {
    int id = -1;
    const QUAD& Q = quad(quad_number);
    for (id=0;id<4;++id) if (Q.line(id)==line_number) break;
    if (id==4) 
      {
	cerr << "  const int TriaContainer<DIM>::find_local_line_index_of_quad(int quad_number, int line_number) const"
	     << endl << "    line " << line_number << " not found in quad " << quad_number << endl;
	abort();
      }
    return id;
    
    
  }
  

  // **************************************************

  template<int DIM>
  const Element<1>& TriaContainer<DIM>::middle_line_of_quad_at_node(int quad_number,int node) const
  {
    const QUAD& Q = quad(quad_number);
    assert(Q.type()==3);

    for (int i=0;i<4;++i)
      {
	const LINE& L = line(quad(Q.child(i)).line((i+1)%4));
	if ((L.node(0)==node)||(L.node(1)==node)) return L;
      }
    abort();
  }

  // **************************************************
  
  template<int DIM>
  const Element<1>& TriaContainer<DIM>::middle_line_of_quad(int quad_number) const
  {
    const QUAD& Q = quad(quad_number);
    assert((Q.type()==1)||(Q.type()==2));
    if (Q.type()==1)      return line(quad(Q.child(0)).line(1));
    else if (Q.type()==2) return line(quad(Q.child(0)).line(2));

    abort();
  }
    
  // **************************************************

  template<int DIM>
  const int TriaContainer<DIM>::line_of_hex(int hn,int ni) const
  {
    assert(DIM==3);
    const HEX& H = hex(hn);

    if (ni<4)
      return line_of_quad_at_nodes(H.quad(0),H.node(ni),H.node((ni+1)%4));

    if ((ni>=4)&&(ni<8))
      return line_of_quad_at_nodes(H.quad(4),H.node(ni),H.node(4+(ni+1)%4));

    if (ni==8)  return line_of_quad_at_nodes(H.quad(1),H.node(1),H.node(5));
    if (ni==9)  return line_of_quad_at_nodes(H.quad(1),H.node(2),H.node(6));
    if (ni==10) return line_of_quad_at_nodes(H.quad(5),H.node(3),H.node(7));
    if (ni==11) return line_of_quad_at_nodes(H.quad(5),H.node(0),H.node(4));
    assert(0);
  }

  // **************************************************

  template<int DIM>
  const int TriaContainer<DIM>::line_of_quad_at_nodes(int qn,int n1,int n2) const
  {
    assert(DIM>=2);
    const QUAD& Q = quad(qn);
    int i=0;
    for (;i<4;++i)
      {
	if (((Q.node(i)==n1)&&(Q.node((i+1)%4)==n2))||
	    ((Q.node(i)==n2)&&(Q.node((i+1)%4)==n1))) break;
      }
    assert(i<4);
    return Q.line(i);
  }
  // **************************************************


  template<int DIM>
  const int TriaContainer<DIM>::child_of_line_at_node(int li,int ni) const
  {
    assert(li>=0);
    const LINE& L = line(li);
    assert(L.nchilds()==2);

    if      (L.node(0)==ni) return L.child(0);
    else if (L.node(1)==ni) return L.child(1);
    else
      {
	cerr << "node " << ni << " not found in line " << li << " nodes:"
	     << line(li).node(0) << "\t" << line(li).node(1) << endl;
	abort();
      }
  }


  // **************************************************

  template<int DIM>
  const int TriaContainer<DIM>::child_of_quad_at_node(int qi,int ni) const
  {
    const QUAD& Q = quad(qi);
    assert(Q.type()>0);

    int node;
    for (node=0;node<4;++node) if (Q.node(node)==ni) break;
    assert(node<4);
    if (Q.type()==3) return Q.child(node);
    if (Q.type()==1) return Q.child(((node+1)%4)/2);
    if (Q.type()==2) return Q.child(node/2);
    
    abort();
  }

  // **************************************************

  template<int DIM>
  void TriaContainer<DIM>::set_adjacent_quad(int li,int id)
  {
    LINE& L = line(li);
    if (L.master()==-1) L.master()=id;
    else if (L.slave()==-1) L.slave()=id;
    else
      {
	cout << "master/slave in line " << li << " both set" << endl;
	cout << L.master() << " " << L.slave() << endl;
	
	assert(0);
      }
  }


  // **************************************************

  template<int DIM>
  void TriaContainer<DIM>::set_adjacent_quad_and_replace(int li,int id, int old)
  {
    LINE& L = line(li);
    
    if (L.master() == old)     L.master() = id;
    else if (L.slave() == old) L.slave() = id;
    else
      {
      if      (L.master()==-1) L.master() = id;
      else if (L.slave()==-1)  L.slave()  = id;
      else 
	{
	  cout << "master/slave in line " << li << " both set and quad " 
	    << old << " not present. quads = " 
	    << L.master() << " " << L.slave() << endl;
	  abort();
	 }
     }
  }


  // **************************************************

  template<int DIM>
  void TriaContainer<DIM>::set_adjacent_hex_and_replace(int li,int _old,int _new)
  {
    assert(DIM==3);
    QUAD& Q = quad(li);
  
    if      (Q.master()==_old) Q.master()=_new;
    else if (Q.slave()==_old)  Q.slave() =_new;
    else if (Q.master()==-1)   Q.master()=_new;
    else if (Q.slave()==-1)    Q.slave() =_new;
    else assert(0);
  }

  // **************************************************
  
  template<int DIM>
  int TriaContainer<DIM>::node_0_of_boundary_quad(const int hex_number,
						  const int boundary_quad) const
  {
    assert(hex_number<nhexes());
    assert(boundary_quad<6);
    
    const HEX&  H = hex(hex_number);
    const QUAD& Q = quad(H.quad(boundary_quad));
    
    int first = 0;
    for (first=0;first<4;++first)
      if (Q.node(0)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad][first]))
	break;
    if (first==4)
      {
	cerr << "node 0 not found in boundary quad " 
	     << boundary_quad << " of hex " << hex_number << endl;
	abort();
      }
    return first;
  }
  
  // **************************************************
  
  template<int DIM>
  int TriaContainer<DIM>::anisotropic_orientation_of_boundary_quad(const int hex_number,const int boundary_quad) const
  {
    assert(hex_number<nhexes());
    assert(boundary_quad<6);

    const HEX&  H = hex(hex_number);
    const QUAD& Q = quad(H.quad(boundary_quad));

    int fn  = node_0_of_boundary_quad(hex_number,boundary_quad);
    
    int dir;
    if (Q.node(1)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad][(fn+1)%4]))
      dir = 0;
    else if (Q.node(1)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad][(fn+3)%4]))
      dir = 1;
    else abort();

    return (fn+dir)%2;
  }
    
  // **************************************************

  template<int DIM>
  pair<int,int> TriaContainer<DIM>::rotation_of_boundary_quad(const int hex_number,
							      const int boundary_quad_number) const
  {
    const HEX&  H = hex(hex_number);
    const QUAD& Q = quad(H.quad(boundary_quad_number));

    // find first node
    int first;
    for (first=0;first<4;++first) 
      if (Q.node(0)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad_number][first])) break;
    assert(first<4);
    // get direction
    int dir=0;
    if (Q.node(1)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad_number][(first+1)%4])) dir=1;
    else if (Q.node(1)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad_number][(first+3)%4])) dir=-1;
    else
      {
	cerr << hex_number << "\t" << boundary_quad_number << endl;
	for (int i=0;i<4;++i) cout << Q.node(i) << " " ; cerr << endl;
	for (int i=0;i<8;++i) cout << H.node(i) << " " ; cerr << endl;
	abort();
      }
    
    // check quad
    if (true)
      {
	if (dir==1)
	  for (int i=0;i<4;++i)
	    assert(Q.node(i)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad_number][(first+i)%4]));
	else if (dir==-1)
	  for (int i=0;i<4;++i)
	    assert(Q.node(i)==H.node(STD_BOUNDARY_QUAD_DIRECTION[boundary_quad_number][(first-i+4)%4]));
      }

    if (dir==-1) first = (first+1)%4;
  
    return make_pair<int,int> (first,dir);
  }



  
  template<> int TriaContainer<1>::refine_cell(const int i,int type) { return refine_line(i,type); }
  template<> int TriaContainer<2>::refine_cell(const int i,int type) { return refine_quad(i,type); }
  template<> int TriaContainer<3>::refine_cell(const int i,int type) { return refine_hex (i,type); }

  template<> int TriaContainer<1>::set_refine_flag(int c, int type) { return set_refine_flag_on_line(c,type); }
  template<> int TriaContainer<2>::set_refine_flag(int c, int type) { return set_refine_flag_on_quad(c,type); }
  template<> int TriaContainer<3>::set_refine_flag(int c, int type) { return set_refine_flag_on_hex (c,type); }

  // ------------------------------
  
  template<int DIM> 
  void TriaContainer<DIM>::clear_refine_flags()
  {
    for (int i=0;i<nlines();++i)            line(i).flag()=0;
    if (DIM>1) for (int i=0;i<nquads();++i) quad(i).flag()=0;
    if (DIM>2) for (int i=0;i<nhexes();++i) hex(i). flag()=0;
  }

  // ------------------------------

  template<int DIM>
  void TriaContainer<DIM>::reinit_hanging_nodes()
  {  
    __hanging_lines.clear();
    __hanging_quads.clear();
    
    if (DIM==2)
      {
	// line hangs, if it's boundary of active quad and refined
	for (int q=0;q<nquads();++q)
	  if (quad(q).type()==0)
	    for (int l=0;l<4;++l)
	      if (line(quad(q).line(l)).type()!=0)
		__hanging_lines[middle_node(line(quad(q).line(l)))] = 
		  quad(q).line(l);
      }
    
    if (DIM==3) // This is isotropic and anisotropic
      {
	// run along all active hexes
	for (int h=0;h<nhexes();++h)
#define H hex(h)
	  if (H.nchilds()==0)
	    {
	      // quads
	      for (int q=0;q<6;++q)
#define Q quad(H.quad(q))
		if (Q.nchilds()>0)
		  {
		    int ty = Q.type();
		    if (ty==3) __hanging_quads[middle_node(Q)]=H.quad(q);
		    else
		      {
			assert(Q.nchilds()==2);
#define Q0 quad(Q.child(0))
#define Q1 quad(Q.child(1))
			if ((Q0.type()==3-ty)&&(Q1.type()==3-ty)) __hanging_quads[middle_node(Q)]=H.quad(q);
#undef Q0
#undef Q1
		      }
		  }
#undef Q
	      
	      // lines
	      for (int l=0;l<12;++l)
		{
		  int li = line_of_hex(h,l);
		  if (line(li).nchilds()>0)
		    {
		      __hanging_lines[middle_node(line(li))]=li;
		    }
		}  
	    }
      }
#undef H
  }

  // ------------------------------

  template<int DIM>
  void TriaContainer<DIM>::reinit_curved()
  {
    typename std::map<int,const BoundaryFunction<DIM>* >::const_iterator it;
    for (it=__curved_boundary.begin();it!=__curved_boundary.end();++it)
      {
	int col                         = it->first;
	const BoundaryFunction<DIM>* BF = it->second;
	HASH_MAP::const_iterator it_b;
	for (it_b = __boundary_lines.begin();it_b!=__boundary_lines.end();++it_b)
	  {
	    if (it_b->second != col) continue;
	    
	    if (line(it_b->first).type()!=0) continue;
	    BF->newton(vertex(line(it_b->first).node(0)));
	    BF->newton(vertex(line(it_b->first).node(1)));
	  }
      }
  }

  // ------------------------------

  template<int DIM>
  void TriaContainer<DIM>::post_refine() 
  {
    reinit_hanging_nodes();
    reinit_curved();
    print_gnuplot("X");

    sanity_test();
  }
  

  // ------------------------------

  template<int DIM> 
  int TriaContainer<DIM>::add(const VERTEX& v)
  { __vertices.push_back(v); return __vertices.size()-1; }


  // ------------------------------

  template<int DIM> 
  int TriaContainer<DIM>::add(const HEX&  h)
  {
    __hexes.push_back(h);
    int id = __hexes.size()-1;
    __hexes[id].id()=id;
    return id;
  }
  
  // ------------------------------

  template<int DIM> 
  int TriaContainer<DIM>::add(const QUAD& q)
  {
    __quads.push_back(q);
    int id = __quads.size()-1;
    __quads[id].id()=id;
    return id;
  }
  
  // ------------------------------

  template<int DIM> 
  int TriaContainer<DIM>::add(const LINE& l)
  {
    __lines.push_back(l);
    int id = __lines.size()-1;
    __lines[id].id()=id;
    return id;
  }

  // ------------------------------

  template<int DIM> 
  void TriaContainer<DIM>::AddBoundaryLine(int l,int color) 
  {
    assert(!line_at_boundary(l));
    __boundary_lines[l]=color;
  }

  // ------------------------------
  
  template<int DIM> 
  void TriaContainer<DIM>::AddBoundaryQuad(int l,int color) 
  {
    assert(!quad_at_boundary(l));
    __boundary_quads[l]=color;
  }
  
  // **************************************************

  template<int DIM> MeshVertex<DIM> TriaContainer<DIM>::new_vertex_on_line(const LINE& L) const
  {
    MeshVertex<DIM> v;
    for (int n=0;n<2;++n) v += vertex(L.node(n));
    v *= 0.5;
    return v;
  }

  // **************************************************

  template<int DIM> MeshVertex<DIM> TriaContainer<DIM>::new_vertex_on_quad(const QUAD& Q) const
  {
    MeshVertex<DIM> v;
    for (int n=0;n<4;++n) v += vertex(Q.node(n));
    v *=0.25;
    return v;
  }

  // **************************************************

  template<int DIM> MeshVertex<DIM> TriaContainer<DIM>::new_vertex_on_hex(const HEX& H) const
  {
    MeshVertex<DIM> v;
    for (int n=0;n<8;++n) v += vertex(H.node(n));
    v *= 0.125;
    return v;
  }

  // **************************************************

  template<int DIM> int TriaContainer<DIM>::set_refine_flag_on_line(int c, int type) 
  {
    int rv=0;
    assert(c<nlines());
    assert(type<2);
    if (type==0) return 0;
    if (!(line(c).flag())) rv=1;
    line(c).flag()|=REF_X;
    return rv;
  }

  // --------------------------------------------------

  template<int DIM> int TriaContainer<DIM>::set_refine_flag_on_quad(int c, int type)
  {
    assert(c<nquads());
    assert(type<4);
#define Q quad(c)
    if (isotropic())
      assert((quad(c).type()==3)||(quad(c).type()==0));
    
    // flag already set?
    if ((Q.flag()&type)==type) return 0;
    if ((Q.type()==type))      return 0;

    Q.flag()|=type;

    if (type&REF_X) 
      {
	set_refine_flag_on_line(Q.line(0),1);
	set_refine_flag_on_line(Q.line(2),1);
      }
    if (type&REF_Y) 
      {
	set_refine_flag_on_line(Q.line(1),1);
	set_refine_flag_on_line(Q.line(3),1);
      }
    
#undef Q

    return 1;
  }

  // --------------------------------------------------
  
  template<int DIM> int TriaContainer<DIM>::set_refine_flag_on_boundary_quad(int h,int bq , int type)
  {
    assert(DIM==3);
    assert(h<nhexes());
    if (type==0) return 0;
    if (type==3) return set_refine_flag_on_quad(hex(h).quad(bq),3);
    
    int nd = anisotropic_orientation_of_boundary_quad(h,bq);
    if (nd==0) return set_refine_flag_on_quad(hex(h).quad(bq),type);
    else return set_refine_flag_on_quad(hex(h).quad(bq),3-type);
  }

  // --------------------------------------------------

  template<int DIM> int TriaContainer<DIM>::set_refine_flag_on_hex(int c, int type)
  {
    assert(c<nhexes());
    if (isotropic())
      assert(type==7);
    if (isotropic())
      assert((hex(c).type()==7)||(hex(c).type()==0));
    
#define H hex(c)

    
    // nothing to do?
    if ((H.flag()&type)==type) return 0;

    H.flag()|=type;


    if (type&REF_X)
      {
	set_refine_flag_on_boundary_quad(c,0,REF_X);
	set_refine_flag_on_boundary_quad(c,2,REF_Y);
	set_refine_flag_on_boundary_quad(c,4,REF_Y);
	set_refine_flag_on_boundary_quad(c,3,REF_Y);
      }
    if (type&REF_Y)
      {
	set_refine_flag_on_boundary_quad(c,0,REF_Y);
	set_refine_flag_on_boundary_quad(c,1,REF_Y);
	set_refine_flag_on_boundary_quad(c,4,REF_X);
	set_refine_flag_on_boundary_quad(c,5,REF_X);
      }
    if (type&REF_Z)
      {
	set_refine_flag_on_boundary_quad(c,1,REF_X);
	set_refine_flag_on_boundary_quad(c,2,REF_X);
	set_refine_flag_on_boundary_quad(c,5,REF_Y);
	set_refine_flag_on_boundary_quad(c,3,REF_X);
      }
    
    return 1;
#undef H  
  }



  // **************************************************


  // **************************************************

  template<int DIM>
  void TriaContainer<DIM>::global_refine()
  {
    clear_refine_flags();

    //CHANGED TO PREREFINE ANISOTROPICALLY

    for (int i=0;i<ncells();++i)
      //      set_refine_flag(i,rand()%3+1);
      set_refine_flag(i);
    refine_cells();
  }

  // **************************************************

  template<int DIM>
  int TriaContainer<DIM>::refine_cells()
  {
    int nref=0;

    first_smooth_flags();
    int ns = 0;
    do
      {
	ns=smooth_flags();
      }
    while(ns>0);

    resolve_flags();

    for (int c=0;c<ncells();++c)
#define C cell(c)
      if (C.type()==0)
	if (C.flag())
	  {
	    if  (!refine_cell(c,C.flag()))
	      abort();
	    C.flag()=0;
	    ++nref;
	  }
#undef C
    clear_refine_flags();

    post_refine();
    
    return nref;
  }


  // **************************************************

  template<int DIM>
  int TriaContainer<DIM>::refine_line(const int line_number, int type)
  {
#define L line(line_number)
    assert(type==1);
    
    if (L.type()==1) return 0;

    VERTEX v = new_vertex_on_line(L);
    int nvi  = add(v);    
    
    // create two new lines
    int l0 = add(LINE(line_number,mesharray<2,int> (-1),getmesharray(L.node(0),nvi)));
    int l1 = add(LINE(line_number,mesharray<2,int> (-1),getmesharray(nvi,L.node(1))));
  
    L.type()=1;
    L.child(0)=l0;
    L.child(1)=l1;

    if (line_at_boundary(line_number)) 
      {
	int color = color_of_boundary_line(line_number);
	AddBoundaryLine(L.child(0),color);
	AddBoundaryLine(L.child(1),color);
      }
    return 1;
#undef L
  } 

  /* -------------------------------------------------- */

  template<int DIM>
  int TriaContainer<DIM>::refine_quad(const int quad_number,int type)
  {
#define Q quad(quad_number)
    if (type==0)        return 0;
    if (Q.type()==type) return 0;
    if (isotropic())
      assert(type==3);
    assert(type<4);

    assert(Q.nchilds()==0);
    
    Q.type()=type;
    
    ////////// REFINE LINES
    if (type&REF_X)
      {
	refine_line(Q.line(0));
	refine_line(Q.line(2));
      }
    if (type&REF_Y)
      {
	refine_line(Q.line(1));
	refine_line(Q.line(3));
      }
    
    
    if (type==REF_X)
      {
	// create new line
	int lid = add(LINE(-1,mesharray<2,int>(-1),getmesharray (middle_node(line(Q.line(0))),middle_node(line(Q.line(2))))));

	
	// create two new quads
	int qid1 = add(QUAD(Q.id(), mesharray<4,int> (-1), getmesharray
			    (Q.node(0), middle_node(line(Q.line(0))), 
			     middle_node(line(Q.line(2))), Q.node(3))));
	int qid2 = add(QUAD(Q.id(), mesharray<4,int> (-1), getmesharray
			    (middle_node(line(Q.line(0))), Q.node(1), 
			     Q.node(2), middle_node(line(Q.line(2))))));

	// set children
	Q.child(0)=qid1;
	Q.child(1)=qid2;
	
	// set lines of childs
	quad(qid1).line(0) = child_of_line_at_node(Q.line(0),Q.node(0));
	quad(qid1).line(1) = lid;
	quad(qid1).line(2) = child_of_line_at_node(Q.line(2),Q.node(3));
	quad(qid1).line(3) = Q.line(3);

	quad(qid2).line(0) = child_of_line_at_node(Q.line(0),Q.node(1));
	quad(qid2).line(1) = Q.line(1);
	quad(qid2).line(2) = child_of_line_at_node(Q.line(2),Q.node(2));
	quad(qid2).line(3) = lid;

	// set master-slaves
	if (DIM==2)
	  {
	    for (int l=0;l<4;++l) set_adjacent_quad_and_replace(quad(qid1).line(l),qid1,quad_number);
	    for (int l=0;l<4;++l) set_adjacent_quad_and_replace(quad(qid2).line(l),qid2,quad_number);
	  }
      }

    if (type==REF_Y)
      {
	// create new line
	int lid = add(LINE(-1,mesharray<2,int>(-1),getmesharray
			   (middle_node(line(Q.line(1))),middle_node(line(Q.line(3))))));
	
	
	// create two new quads
	int qid1 = add(QUAD(Q.id(), mesharray<4,int> (-1), getmesharray
			    (Q.node(0), Q.node(1), middle_node(line(Q.line(1))), 
			     middle_node(line(Q.line(3))))));
	int qid2 = add(QUAD(Q.id(), mesharray<4,int> (-1), getmesharray
			    (middle_node(line(Q.line(3))), 
			     middle_node(line(Q.line(1))), Q.node(2), Q.node(3))) );
	
	// set children
	Q.child(0)=qid1;
	Q.child(1)=qid2;
	
	// set lines of childs
	quad(qid1).line(0) = Q.line(0);
	quad(qid1).line(1) = child_of_line_at_node(Q.line(1),Q.node(1));
	quad(qid1).line(2) = lid;
	quad(qid1).line(3) = child_of_line_at_node(Q.line(3),Q.node(0));


	quad(qid2).line(0) = lid;
	quad(qid2).line(1) = child_of_line_at_node(Q.line(1),Q.node(2));
	quad(qid2).line(2) = Q.line(2);
	quad(qid2).line(3) = child_of_line_at_node(Q.line(3),Q.node(3));

	// set master-slaves
	if (DIM==2)
	  {
	    for (int l=0;l<4;++l) set_adjacent_quad_and_replace(quad(qid1).line(l),qid1,quad_number);
	    for (int l=0;l<4;++l) set_adjacent_quad_and_replace(quad(qid2).line(l),qid2,quad_number);
	  }
      }
    
    // create new middle vertex	  
    if (type==(REF_X|REF_Y))
      {
	VERTEX v = new_vertex_on_quad(Q);
	int nvi  = add(v);
    
	// create 3 lines
	/**
	 *   -------
	 *  |   2   |
	 *  |-3-|-1-|
	 *  |   0   |
	 *   -------
	 **/
	// id's of new lines
	mesharray<4,int> lid;
	for (int i=0;i<4;++i)
	  lid[i]=add(LINE(-1,mesharray<2,int>(-1),getmesharray(middle_node(line(Q.line(i))),nvi)));
	
	// ***** 4)
	// create the 4 quads
	// id of new quads
	mesharray<4,int> cid;
	for (int i=0;i<4;++i)
	  {
	    // vertices
	    mesharray<4,int> vv;
	    vv[i]=Q.node(i);
	    vv[(i+1)%4]=middle_node(line(Q.line(i)));
	    vv[(i+2)%4]=nvi;
	    vv[(i+3)%4]=middle_node(line(Q.line((i+3)%4)));
	    
	    cid[i]=add(QUAD(Q.id(),mesharray<4,int>(-1),vv));
	    	    
	    // ***** 5)
	    // set lines of childs
	    quad(cid[i]).line(i)=child_of_line_at_node(Q.line(i),vv[i]);
	    quad(cid[i]).line((i+1)%4)=lid[i];
	    quad(cid[i]).line((i+2)%4)=lid[(i+3)%4];
	    quad(cid[i]).line((i+3)%4)=child_of_line_at_node(Q.line((i+3)%4),vv[i]);
	    // ***** 6)
	    // set master/slave-info of lines
	    if (DIM==2)
	      {
		set_adjacent_quad(child_of_line_at_node(Q.line(i),vv[i]),cid[i]);
		set_adjacent_quad(lid[i],cid[i]);
		set_adjacent_quad(lid[(i+3)%4],cid[i]);
		set_adjacent_quad(child_of_line_at_node(Q.line((i+3)%4),vv[i]),cid[i]);
	      }
	    Q.child(i)=cid[i];
	  }
      }
    
    // set type of new quads to zero
    for (int c=0;c<Q.nchilds();++c)
      quad(Q.child(c)).type()=0;
    
    // set boundary quads
    if (DIM==3)
      if (quad_at_boundary(quad_number))
	{
	  int color = color_of_boundary_quad(quad_number);
	  for (int c=0;c<Q.nchilds();++c)
	    AddBoundaryQuad(Q.child(c),color);
	}
#undef Q
    return 1;
  }

  template<int DIM>
  void TriaContainer<DIM>::set_boundary_quads(int hn,int q0,int q1,int q2,int q3,int q4,int q5)
  {
#define H hex(hn)
    H.quad(0)=q0;
    H.quad(1)=q1;
    H.quad(2)=q2;
    H.quad(3)=q3;
    H.quad(4)=q4;
    H.quad(5)=q5;
    
    // check them
    for (int i=0;i<6;++i)
      node_0_of_boundary_quad(hn,i);
#undef H
  }
  

  template<int DIM>
  int TriaContainer<DIM>::split_hex_isotropic(const int hn, const mesharray<27,int>& cn, int type)
  {
#define H hex(hn)
#define Q0 quad(hex(hn).quad(0))
#define Q1 quad(hex(hn).quad(1))
#define Q2 quad(hex(hn).quad(2))
#define Q3 quad(hex(hn).quad(3))
#define Q4 quad(hex(hn).quad(4))
#define Q5 quad(hex(hn).quad(5))
    
    assert(type==7);

    // 6 new lines
    LINE L;
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[4], cn[13])); int li0 = add (L);
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[14],cn[13])); int li1 = add (L);
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[16],cn[13])); int li2 = add (L);
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[10],cn[13])); int li3 = add (L);
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[22],cn[13])); int li4 = add (L);
    L.init(-1,mesharray<2,int> (-1), getmesharray(cn[12],cn[13])); int li5 = add (L);
    
    // 12 new quads
    QUAD Q;
    Q.new_by_lines(line(quad(Q.child(0)).line(1)),
		   line(quad(Q.child(0)).line(1)),
		   line(quad(Q.child(0)).line(1)),
		   line(quad(Q.child(0)).line(1)));
    
    
#undef H
#undef Q0
#undef Q1
#undef Q2
#undef Q3
#undef Q4
#undef Q5
  }
  


  template<int DIM>
  int TriaContainer<DIM>::split_hex_one_way(const int hn, const mesharray<27,int>& cn, int type)
  {
    assert(hn<nhexes());
    assert((type==REF_X)||(type==REF_Y)||(type==REF_Z));
#define H hex(hn)
    assert((H.type()==REF_X)||(H.type()==REF_Y)||(H.type()==REF_Z));
    QUAD Q;
    
#define mloq middle_line_of_quad
    if (type==REF_X)      Q.new_by_lines(mloq(H.quad(0)),mloq(H.quad(2)),mloq(H.quad(4)),mloq(H.quad(3)));
    else if (type==REF_Y) Q.new_by_lines(mloq(H.quad(0)),mloq(H.quad(1)),mloq(H.quad(4)),mloq(H.quad(5)));
    else if (type==REF_Z) Q.new_by_lines(mloq(H.quad(1)),mloq(H.quad(2)),mloq(H.quad(5)),mloq(H.quad(3)));
    else abort();
#undef mloq

    int nq = add(Q);

    int h1=-1,h2=-1;

    if (type==REF_X)
      {
	h1 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[0],cn[1],cn[7],cn[6],cn[18],cn[19],cn[25],cn[24])));
	h2 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[1],cn[2],cn[8],cn[7],cn[19],cn[20],cn[26],cn[25])));
      }
    if (type==REF_Y)
      {
	h1 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[0],cn[2],cn[5],cn[3],cn[18],cn[20],cn[23],cn[21])));
	h2 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[3],cn[5],cn[8],cn[6],cn[21],cn[23],cn[26],cn[24])));
      }
    if (type==REF_Z)
      {
	h1 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[0],cn[2],cn[8],cn[6],cn[9],cn[11],cn[17],cn[15])));
	h2 = add(HEX(hn, mesharray<8,int> (-1), getmesharray
		     (cn[9],cn[11],cn[17],cn[15],cn[18],cn[20],cn[26],cn[24])));
      }
    H.child(0)=h1;
    H.child(1)=h2;
    
    // boundaries
#define coqan child_of_quad_at_node
    if (type==REF_X)
      {
	set_boundary_quads(h1, coqan(H.quad(0),cn[0]),nq,coqan(H.quad(2),cn[6]),
			   coqan(H.quad(3),cn[0]),coqan(H.quad(4),cn[18]),H.quad(5));
	set_boundary_quads(h2, coqan(H.quad(0),cn[2]),H.quad(1),coqan(H.quad(2),cn[8]),
			   coqan(H.quad(3),cn[2]),coqan(H.quad(4),cn[20]),nq);
      }
    if (type==REF_Y)
      {
	set_boundary_quads(h1, coqan(H.quad(0),cn[0]),coqan(H.quad(1),cn[2]),nq,
			   H.quad(3),coqan(H.quad(4),cn[18]),coqan(H.quad(5),cn[0]));
	set_boundary_quads(h2, coqan(H.quad(0),cn[6]),coqan(H.quad(1),cn[8]),H.quad(2),
			   nq,coqan(H.quad(4),cn[24]),coqan(H.quad(5),cn[6]));
      }
    if (type==REF_Z)
      {
	set_boundary_quads(h1, H.quad(0),coqan(H.quad(1),cn[2]),coqan(H.quad(2),cn[6]),
			   coqan(H.quad(3),cn[0]),nq,coqan(H.quad(5),cn[0]));
	set_boundary_quads(h2, nq,coqan(H.quad(1),cn[20]),coqan(H.quad(2),cn[24]),
			   coqan(H.quad(3),cn[20]),H.quad(4),coqan(H.quad(5),cn[18]));
      }
#undef coqan

    // adjacency
    for (int q=0;q<6;++q)
      set_adjacent_hex_and_replace(hex(h1).quad(q),hn,h1);
    for (int q=0;q<6;++q)
      set_adjacent_hex_and_replace(hex(h2).quad(q),hn,h2);

#undef H
    return 1;
  }
  
  
  template<int DIM>
  int TriaContainer<DIM>::refine_hex(const int hex_number,int type)
  {
#define H hex(hex_number)
    assert(type<8);
    
    if (isotropic())
      assert(type==7);
    
    // nothing to do?
    if (type==0)
      return 0;
    
    // already refined in that direction
    if (type==H.type())
      return 0;
    
    // No re refinement possible
    if (H.type()>0) abort();

    // set type of hex
    H.type()=type;


    /////// REFINE BOUNDARIES
    if (type==(REF_X|REF_Y|REF_Z))
      {
	for (int i=0;i<6;++i)
	  refine_boundary_quad(hex_number,i,3);
      }
    else assert(0);



    
    //////////////////////////////
    ////////////////////////////// GET INDICES OF ALL NODES
    //////////////////////////////

    mesharray<27,int> cn(-1);
    // nodes in corners
    cn[0] = H.node(0);    cn[2] = H.node(1);
    cn[6] = H.node(3);    cn[8] = H.node(2);
    cn[18]= H.node(4);    cn[20]= H.node(5);
    cn[24]= H.node(7);    cn[26]= H.node(6);

    // nodes on lines
    if (type&REF_X)
      {
	cn[1] =middle_node(line(line_of_hex(hex_number,0)));
	cn[7] =middle_node(line(line_of_hex(hex_number,2)));
	cn[19]=middle_node(line(line_of_hex(hex_number,4)));
	cn[25]=middle_node(line(line_of_hex(hex_number,6)));
      }
    if (type&REF_Y)
      {
	cn[3] =middle_node(line(line_of_hex(hex_number,3)));
	cn[5] =middle_node(line(line_of_hex(hex_number,1)));
	cn[21]=middle_node(line(line_of_hex(hex_number,7)));
	cn[23]=middle_node(line(line_of_hex(hex_number,5)));
      }
    if (type&REF_Z)
      {
	cn[9] =middle_node(line(line_of_hex(hex_number,11)));
	cn[11]=middle_node(line(line_of_hex(hex_number,8)));
	cn[15]=middle_node(line(line_of_hex(hex_number,10)));
	cn[17]=middle_node(line(line_of_hex(hex_number,9)));
      }
    // nodes on quads
    if ((type&(REF_X|REF_Y)) == (REF_X|REF_Y))
      {
	cn[4] = middle_node(quad(H.quad(0)));
	cn[22]= middle_node(quad(H.quad(4)));
      }
    if ((type&(REF_X|REF_Z)) == (REF_X|REF_Z))
      {    
	cn[10]= middle_node(quad(H.quad(3)));
	cn[16]= middle_node(quad(H.quad(2)));
      }
    if ((type&(REF_Y|REF_Z)) == (REF_Y|REF_Z))
      {    
	cn[12]= middle_node(quad(H.quad(5)));
	cn[14]= middle_node(quad(H.quad(1)));
      }
    // ISOTROPIC REFINEMENT: create new middle node
    if (type==(REF_X|REF_Y|REF_Z))
      {
	VERTEX v;
	for (int i=0;i<8;++i) v+=vertex(H.node(i));
	v*=0.125;
	cn[13] = add(v);
      }
    


    //////////////////////////////
    //////////////////////////////  DIFFERENT CASES
    //////////////////////////////

    // isotropic: create 6 new line, 12 new quads, 6 hexes
    if (type==(REF_X|REF_Y|REF_Z))
      split_hex_isotropic(hex_number,cn,type);
    
      {
	// create new lines (ordered with quads)

	// QUAD Q;
	// // front quads
	// Q.new_by_lines(mloqan(H.quad(0),cn[1]), line(li0),line(li3),mloqan(H.quad(3),cn[1]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(0),cn[5]), line(li0),line(li1),mloqan(H.quad(1),cn[5]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(0),cn[7]), line(li0),line(li2),mloqan(H.quad(2),cn[7]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(0),cn[3]), line(li0),line(li5),mloqan(H.quad(5),cn[3]));
	// new_quads.push_back(add(Q));	
	// // back quads
	// Q.new_by_lines(mloqan(H.quad(4),cn[19]), line(li4),line(li3),mloqan(H.quad(3),cn[19]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(4),cn[23]), line(li4),line(li1),mloqan(H.quad(1),cn[23]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(4),cn[25]), line(li4),line(li2),mloqan(H.quad(2),cn[25]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(4),cn[21]), line(li4),line(li5),mloqan(H.quad(5),cn[21]));
	// new_quads.push_back(add(Q));	
	// // side quads
	// Q.new_by_lines(mloqan(H.quad(1),cn[11]), line(li1),line(li3),mloqan(H.quad(3),cn[11]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(2),cn[17]), line(li2),line(li1),mloqan(H.quad(1),cn[17]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(5),cn[15]), line(li5),line(li2),mloqan(H.quad(2),cn[15]));
	// new_quads.push_back(add(Q));
	
	// Q.new_by_lines(mloqan(H.quad(3),cn[9]),  line(li3),line(li5),mloqan(H.quad(5),cn[9]));
	// new_quads.push_back(add(Q));	
    
      }
    

    ////////// The isotropic case. New lines in the middle are necessary

    
    if ((type==REF_X)||(type==REF_Y)||(type==REF_Z))
      split_hex_one_way(hex_number,cn,type);

    else abort();
    
    
//     if (type==REF_X)
//       {
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(0,1,7,6,18,19,25,24),
// 		       getmesharray(H.quad(0),nq,coqan(H.quad(2),cn[6])),coqan(
// 		       );
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(1,2,8,7,19,20,26,25),
// 		       );
//       }
//     if (type==REF_Y)
//       {
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(0,2,5,3,18,20,23,21),
// 		       );
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(3,5,8,6,21,23,26,24),
// 		       );
//       }
//     if (type==REF_Z)
//       {
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(0,2,8,6,9,11,17,15),
// 		       );
// 	create_new_hex(hex_number,cn,
// 		       getmesharray(9,11,17,15,18,20,26,24),
// 		       );
//       }
    

// #define mloq middle_line_of_quad
// #define mloqan middle_line_of_quad_at_node

//     // create new lines (ordered with quads)
//     LINE L;
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[4], cn[13])); int li0 = add (L);
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[14],cn[13])); int li1 = add (L);
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[16],cn[13])); int li2 = add (L);
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[10],cn[13])); int li3 = add (L);
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[22],cn[13])); int li4 = add (L);
//     L.init(-1,mesharray<2,int> (-1), getmesharray(cn[12],cn[13])); int li5 = add (L);

//     // create new quads in the middle of the cell
//     vector<int> new_quads;
    
    
//     // new quads (ordered with lines)

  
// #undef mloq
// #undef mloqan

// #define coqan child_of_quad_at_node
// #define coqal child_of_quad_at_line
//     // create new children
  
//     assert(new_quads.size()==12);
      
//     int tmp = create_child_hex(hex_number,cn,getmesharray(0,1,4,3,9,10,13,12)     ,
// 			       getmesharray(coqan(H.quad(0),cn[0]),
// 					   new_quads[0]             ,
// 					   new_quads[3]             ,
// 					   coqan(H.quad(3),cn[0]),new_quads[11],
// 					   coqan(H.quad(5),cn[0]) ));
//     H.child(0) = tmp;
//     tmp = create_child_hex(hex_number,cn,getmesharray(1,2,5,4,10,11,14,13)    ,
// 			   getmesharray(coqan(H.quad(0),cn[2]),
// 				       coqan(H.quad(1),cn[2]),
// 				       new_quads[1]             ,
// 				       coqan(H.quad(3),cn[2]),
// 				       new_quads[8], 
// 				       new_quads[0]));
//     H.child(1) = tmp;
//     tmp= create_child_hex(hex_number,cn,getmesharray(4,5,8,7,13,14,17,16)    ,
// 			  getmesharray(coqan(H.quad(0),cn[8]),
// 				      coqan(H.quad(1),cn[8]),
// 				      coqan(H.quad(2),cn[8]),
// 				      new_quads[1]             ,
// 				      new_quads[9], 
// 				      new_quads[2]));
//     H.child(2) = tmp;
//     tmp = create_child_hex(hex_number,cn,getmesharray(3,4,7,6,12,13,16,15)    ,
// 			   getmesharray(coqan(H.quad(0),cn[6]),
// 				       new_quads[2]             ,
// 				       coqan(H.quad(2),cn[6]),
// 				       new_quads[3]             ,
// 				       new_quads[10],
// 				       coqan(H.quad(5),cn[6]) ));
//     H.child(3) = tmp;
      
//     tmp = create_child_hex(hex_number,cn,getmesharray(9,10,13,12,18,19,22,21) ,
// 			   getmesharray(new_quads[11],
// 				       new_quads[4]              ,
// 				       new_quads[7]              ,
// 				       coqan(H.quad(3),cn[18]),
// 				       coqan(H.quad(4),cn[18]),
// 				       coqan(H.quad(5),cn[18]) ));
//     H.child(4) = tmp;
//     tmp = create_child_hex(hex_number,cn,getmesharray(10,11,14,13,19,20,23,22),
// 			   getmesharray(new_quads[ 8],
// 				       coqan(H.quad(1),cn[20]),
// 				       new_quads[5]              ,
// 				       coqan(H.quad(3),cn[20]),
// 				       coqan(H.quad(4),cn[20]),
// 				       new_quads[4] ));
//     H.child(5) = tmp;
//     tmp = create_child_hex(hex_number,cn,getmesharray(13,14,17,16,22,23,26,25),
// 			   getmesharray(new_quads[ 9],
// 				       coqan(H.quad(1),cn[26]),
// 				       coqan(H.quad(2),cn[26]),
// 				       new_quads[5]             ,
// 				       coqan(H.quad(4),cn[26]),
// 				       new_quads[6] ));
//     H.child(6) = tmp;
      
//     tmp = create_child_hex(hex_number,cn,getmesharray(12,13,16,15,21,22,25,24),
// 			   getmesharray(new_quads[10],
// 				       new_quads[6]              ,
// 				       coqan(H.quad(2),cn[24]),
// 				       new_quads[7]             ,
// 				       coqan(H.quad(4),cn[24]),
// 				       coqan(H.quad(5),cn[24]) ));
//     H.child(7) = tmp;

// #undef coqan    
// #undef coqal
  
  
    return 1;
  
#undef H
  }

  
  
  template<int DIM>
  int TriaContainer<DIM>::refine_boundary_quad(int hn, int qn, int type)
  {
    assert(DIM==3);
    
    if ((type==0)||(type==3)) return refine_quad(hex(hn).quad(qn),type);
    
    int nd = anisotropic_orientation_of_boundary_quad(hn,qn);
    if (nd==0) return refine_quad(hex(hn).quad(qn),type);
    else return       refine_quad(hex(hn).quad(qn),3-type);
    
  }
  

  // **************************************************

  template<int DIM>
  int TriaContainer<DIM>::create_child_hex(const int hex_number,const mesharray<27,int>& cnode,
					   const mesharray<8,int>& n,
					   const mesharray<6,int>& q)
  {
    for (int i=0;i<8;++i) 
      {
	if (cnode[n[i]]==-1) cout << i << "\t" << n[i] << endl;
	assert(cnode[n[i]]!=-1);
      }
    for (int i=0;i<8;++i) assert(cnode[n[i]]<nvertices());
    for (int i=0;i<6;++i) assert(q[i]!=-1);
    for (int i=0;i<6;++i) assert(q[i]<nquads());
  

    HEX H;
  
    mesharray<8,int> childs(-1),nodes;
    // set nodes
    for (int i=0;i<8;++i)  nodes[i]=cnode[n[i]]; 
    // set quads
    for (int i=0;i<6;++i)  H.quad(i) = q[i];
  
    H.init(hex_number,childs, nodes);
    int hn = add(H);
    // set adjacent info in quads (replace number if father is present)
    for (int i=0;i<6;++i)
      set_adjacent_hex_and_replace(H.quad(i),hex_number,hn);
  
    return hn;
  }

  

  void SANITY_EXIT(const string& s)
  { cerr << "Sanity: " << s << endl; abort(); }

  void SANITY_EXIT(const string& s, const int& i) 
  { cerr << "Sanity: " << s << " " << i << endl; abort(); }

  void SANITY_EXIT(const string& s, const int& i1, const int& i2) 
  { cerr << "Sanity: " << s << " " << i1 << " " << i2 << endl; abort(); }



  // --------------------------------------------------

  template<int DIM>
  void TriaContainer<DIM>::sanity_test() const
  {
    if (DIM>=2)
      sanity_test_quads();
    if (DIM==3)
      sanity_test_hexes();
    sanity_test_hanging_nodes();
  }


  // --------------------------------------------------

  template<int DIM>
  void TriaContainer<DIM>::sanity_test_lines() const{}


  // --------------------------------------------------

  template<int DIM>
  void TriaContainer<DIM>::sanity_test_quads() const
  {
    if (isotropic())
      for (int h=0;h<nquads();++h)
	assert((quad(h).type()==0)||(quad(h).type()==3));
  }



  // --------------------------------------------------

  template<int DIM>
  void TriaContainer<DIM>::sanity_test_hexes() const
  {
    if (isotropic())
      for (int h=0;h<nhexes();++h)
	assert((hex(h).type()==0)||(hex(h).type()==7));
  }


  // --------------------------------------------------

  template<int DIM>
  void TriaContainer<DIM>::sanity_test_hanging_nodes() const
  {
    HASH_MAP::const_iterator hl;
    HASH_MAP::const_iterator hq;

    for (hl = __hanging_lines.begin();hl!=__hanging_lines.end();++hl)
      for (int n=0;n<2;++n)
	if (node_is_hanging(line(hl->second).node(n)))
	  SANITY_EXIT("Hanging node # hangs on hanging node.", line(hl->second).node(n));

    if (DIM==3)  
      for (hq = __hanging_quads.begin();hq!=__hanging_quads.end();++hq)
	for (int n=0;n<4;++n)
	  if (node_is_hanging(quad(hq->second).node(n)))
	    SANITY_EXIT("Hanging node # hangs on hanging node.", quad(hq->second).node(n));
  }


  


  template<int DIM>
  void TriaContainer<DIM>::print_gnuplot(string fname,int extra) const
  {
    
    stringstream data;
    data << fname << "_data";
    
    ofstream OUT(data.str().c_str());
    for (int i=0;i<nlines();++i) 
      if (line(i).type()==0)
	OUT << vertex(line(i).node(0)) << endl
	    << vertex(line(i).node(1)) << endl << endl << endl;
    OUT.close();

    double minx=0,maxx=0,miny=0,maxy=0,maxz=0,minz=0;
    minx=maxx=vertex(0).x();
    miny=maxy=vertex(0).y();
    if (DIM==3)
      minz=maxz=vertex(0).z();
    
    for (int i=0;i<nvertices();++i)
      {
	minx = std::min(minx, vertex(i).x());
	maxx = std::max(maxx, vertex(i).x());
	miny = std::min(miny, vertex(i).y());
	maxy = std::max(maxy, vertex(i).y());
	if (DIM==3)
	  {
	    minz = std::min(minz, vertex(i).z());
	    maxz = std::max(maxz, vertex(i).z());
	  }
      }
    double dx = maxx-minx;    assert(dx>0);
    double dy = maxx-minx;    assert(dy>0);
    double dz = maxz-minz;    if (DIM==3) assert(dz>0);
    
    OUT.open(fname.c_str());
    OUT << "set xrange [" << minx-0.05*dx << ":" << maxx+0.05*dx << "]" << endl;
    OUT << "set yrange [" << miny-0.05*dy << ":" << maxy+0.05*dy << "]" << endl;
    if (DIM==3)
      OUT << "set zrange [" << minz-0.05*dz << ":" << maxz+0.05*dz << "]" << endl;
    OUT << "set sty da li" << endl;
    

    // print cell-id
    if (extra&1)
      {
	for (int c=0;c<ncells();++c)
	  //	  if (cell(c).type()==0)
	    {
	      MeshVertex<DIM> v;
	      for (int k=0;k<4;++k) 
		v += vertex(cell(c).node(k));
	      v *= 0.25;

	      OUT << "set label \"" << c << "\" at " << v.x() << "," << v.y() << endl;
	    }
      }

    // print line-id
    if (extra&2)
      {
	for (int c=0;c<nlines();++c)
	  //	  if (line(c).type()==0)
	    {
	      MeshVertex<DIM> v;
	      for (int k=0;k<2;++k) 
		v += vertex(line(c).node(k));
	      v *= 0.5;

	      OUT << "set label \"" << c << "\" at " << v.x() << "," << v.y() << endl;
	    }
      }
    

    if (DIM==2)
      OUT << "plot \"" << data.str() << "\"" << endl;
    else
      OUT << "splot \"" << data.str() << "\"" << endl;
    
    OUT << "pause -1" << endl;
    OUT.close();
  }
  

  // **************************************************

  template<>
  void TriaContainer<2>::read_inp(const string& fname)
  {
    reset();
  
    ifstream IN(fname.c_str());
    if (!IN)
      {
	cerr << "could'nt open mesh-file " << fname << endl;
	assert(0);
      }
    int nnodes,nothers,d1,d2,d3,index;
    double dv;
  
    IN >> nnodes >> nothers >> d1 >>d2 >> d3;
    assert(IN);
    assert((d1==0)&&(d2==0)&&(d3==0));

    // read nodes
    MeshVertex<2> v;
    for (int i=0;i<nnodes;++i)
      {
	IN >> index >> v >> dv;
	assert(IN);
	__vertices.push_back(v);
      }

    // read quads & lines
    string str;
    int color;
    mesharray<4,int> quad_nodes;
    mesharray<2,int> line_nodes;
    for (int i=0;i<nothers;++i)
      {
	IN >> index >> color >> str;
	// at the moment no boundaries
	if (str=="quad")
	  {
	    IN >> quad_nodes;
	    QUAD Q;
	    Q.init(-1,mesharray<4,int> (-1),quad_nodes);
	    Q.type()=0;
	    add(Q);
	  
	    assert(IN);
	  }
	else if (str=="line")
	  {
	    IN >> line_nodes;
	    LINE L;
	    L.init(-1,mesharray<2,int> (-1),line_nodes);
	    L.type()=0;
	    int ni = add(L);
	    AddBoundaryLine(ni,color);
	  }
      }

    // create inner lines
    for (int i=0;i<nquads();++i)
      {
	QUAD& Q = quad(i);
      
	for (int j=0;j<4;++j)
	  {
	    mesharray<2,int> n;
	    n[0]=Q.node(j);
	    n[1]=Q.node((j+1)%4);
	    LINE L;
	    L.init(-1,mesharray<2,int> (-1),n);
	    int li = find(L);
	    // line not present
	    if (li==-1)
	      {
		int ll = add(L);
	      
		Q.line(j)=ll;
		line(ll).master()=Q.id();
	      }
	    // line already present
	    else
	      {
		if (line(li).master()==-1) line(li).master()=Q.id();
		else if (line(li).slave()==-1) line(li).slave()=Q.id();
		else assert(0);
		Q.line(j)=li;
	      }
	  }
      
      }
  }

  // **************************************************

  template<>
  void TriaContainer<3>::read_inp(const string& fname)
  {
    reset();
  
    ifstream IN(fname.c_str());
    if (!IN)
      {
	cerr << "could'nt open mesh-file " << fname << endl;
	assert(0);
      }
    int nnodes,nothers,d1,d2,d3,index;

    IN >> nnodes >> nothers >> d1 >>d2 >> d3;
    assert(IN);
    assert((d1==0)&&(d2==0)&&(d3==0));

    // read nodes
    MeshVertex<3> v;
    for (int i=0;i<nnodes;++i)
      {
	IN >> index >> v ;
	assert(IN);
	__vertices.push_back(v);
      }

    // read quads & lines
    string str;
    int color;
    mesharray<8,int> nodes;
    mesharray<4,int> bnodes;
    for (int i=0;i<nothers;++i)
      {
	IN >> index >> color >> str;
	// at the moment no boundaries

	if (str=="hex")
	  {
	    IN >> nodes;
	    HEX H;
	    H.init(-1,mesharray<8,int> (-1),nodes);
	    H.type()=0;
	    add(H);
	    assert(IN);
	  }
	else if (str=="quad")
	  {
	    IN >> bnodes;

	    // create four boundary lines
	    mesharray<4,int> bli;
	    for (int bl=0;bl<4;++bl)
	      {
		mesharray<2,int> bln;
		bln[0]=bnodes[bl];
		bln[1]=bnodes[(bl+1)%4];
		LINE BL;
		BL.init(-1,mesharray<2,int> (-1),bln);
		bli[bl] = find(BL);
		if (bli[bl]==-1) 
		  {
		    bli[bl] = add(BL);
		    AddBoundaryLine(bli[bl],color);
		  }
	      }
	  

	    QUAD BQ;
	    BQ.init(-1,mesharray<4,int> (-1),bnodes);
	    BQ.type()=0;
	    for (int bl=0;bl<4;++bl)
	      BQ.line(bl)= bli[bl];
	    
	    int ni = add(BQ);
	    AddBoundaryQuad(ni,color);
	    assert(IN);
	  }
      
      }

    // create inner quads
    for (int i=0;i<nhexes();++i)
      {
	HEX& H = hex(i);

	int fid[6][4]= {{0,1,2,3}, {1,5,6,2}, {3,2,6,7},
			{0,1,5,4}, {4,5,6,7}, {0,4,7,3}};
      
	for (int j=0;j<6;++j)
	  {
	    mesharray<4,int> n;
	    for (int jj=0;jj<4;++jj) n[jj]=H.node(fid[j][jj]);
	    QUAD Q;
	    Q.init(-1,mesharray<4,int> (-1),n);
	    int li = find(Q);
	    // line not present
	    if (li==-1)
	      {
		int ll = add(Q);
	      
		H.quad(j)=ll;
		quad(ll).master()=H.id();
	      }
	    // line already present
	    else
	      {
		if (quad(li).master()==-1) quad(li).master()=H.id();
		else if (quad(li).slave()==-1) quad(li).slave()=H.id();
		else assert(0);
		H.quad(j)=li;
	      }
	  }
      }
    // create inner lines, without master/slave info
    for (int i=0;i<nquads();++i)
      {
	QUAD& Q = quad(i);
	for (int j=0;j<4;++j)
	  {
	    mesharray<2,int> n;
	    n[0]=Q.node(j);
	    n[1]=Q.node((j+1)%4);
	    LINE L;
	    L.init(-1,mesharray<2,int> (-1),n);
	    int li = find(L);
	    // line not present
	    if (li==-1)
	      {
		int ll = add(L);
	      
		Q.line(j)=ll;
	      }
	    // line already present
	    else
	      {
		Q.line(j)=li;
	      }
	  }
      }  
  }


  template<int DIM>
  void TriaContainer<DIM>::write_gup(const std::string& _fn) const
  {
    string fn = _fn;
    int name_size = fn.size();
    if(name_size<4) fn += ".gup";
    else if(fn.substr(name_size-4,4)!=".gup") fn += ".gup";
  
    ofstream OUT(fn.c_str());
  
    if (!OUT)
      {
	cerr << "TriaContainer::write_gup: could not open file " << fn << " for writing" << endl;
	abort();
      }

    // verteces
    OUT << __vertices.size() << endl;
    for (int i=0;i<__vertices.size();++i) OUT << __vertices[i] << endl;
    // hexes
    OUT << __hexes.size() << endl;
    for (int i=0;i<__hexes.size();++i) OUT << __hexes[i] << endl;
    // quads
    OUT << __quads.size() << endl;
    for (int i=0;i<__quads.size();++i) OUT << __quads[i] << endl;
    // lines
    OUT << __lines.size() << endl;
    for (int i=0;i<__lines.size();++i) OUT << __lines[i] << endl;

    // hanging & boundary
    this->write_data(__hanging_lines,OUT);
    this->write_data(__hanging_quads,OUT);
    this->write_data(__boundary_lines,OUT);
    this->write_data(__boundary_quads,OUT);
  
    OUT.close();
    //   cerr << "write mesh:" << endl
    //        << "\t" << __vertices.size() << " " << __hexes.size() << " " << __quads.size() << " " 
    //        << __lines.size()<< endl;
    //   cerr << "\t" << __hanging_lines.size() << " "
    //        << __hanging_quads.size() << " "
    //        << __boundary_lines.size() << " "
    //        << __boundary_quads.size() << endl;

  }


  template<int DIM>
  void TriaContainer<DIM>::read_gup(const std::string& _fn) 
  {
    string fn = _fn;
    int name_size = fn.size();
    if(name_size<4) fn += ".gup";
    else if(fn.substr(name_size-4,4)!=".gup") fn += ".gup";

    ifstream IN(fn.c_str());
    if (!IN)
      {
	cerr << "TriaContainer::read_gup: could not open file " << fn << " for reading" << endl;
	abort();
      }
    reset();


    int tmp;
    VERTEX tmpv;
    // vertices
    IN >> tmp; if (tmp>0) __vertices.resize(tmp);
    for (int i=0;i<__vertices.size();++i) IN >> __vertices[i];
    // heses
    IN >> tmp; if (tmp>0) __hexes.resize(tmp);
    for (int i=0;i<__hexes.size();++i) IN >> __hexes[i];
    // quads
    IN >> tmp; if (tmp>0) __quads.resize(tmp);
    for (int i=0;i<__quads.size();++i) IN >> __quads[i];
    // lines
    IN >> tmp; if (tmp>0) __lines.resize(tmp);
    for (int i=0;i<__lines.size();++i) IN >> __lines[i];
  

    // hanging & boundary
    this->read_data(__hanging_lines,IN);
    this->read_data(__hanging_quads,IN);
    this->read_data(__boundary_lines,IN);
    this->read_data(__boundary_quads,IN);
  
    assert(IN);
  
    IN.close();
  }

  //////////////////////////////////////////////////
  ////////////////////////////////////////////////// binary
  //////////////////////////////////////////////////


  template<int DIM>
  void TriaContainer<DIM>::write_gip(const std::string& _fn) const
  {
    string fn = _fn;
    int name_size = fn.size();
    if(name_size<4) fn += ".gip";
    else if(fn.substr(name_size-4,4)!=".gip") fn += ".gip";
  
    ofstream OUT(fn.c_str(),ios::binary|ios::out);
  
    if (!OUT)
      {
	cerr << "TriaContainer::write_gip: could not open file " << fn << " for writing" << endl;
	abort();
      }

    // verteces
    OUT << __vertices.size() << endl;

    OUT << "[";
    OUT.write(reinterpret_cast<const char*>(&__vertices[0][0]), sizeof(double) * DIM * __vertices.size());
    OUT << "]" << endl;

    // hexes
    OUT << __hexes.size() << endl;
    for (int i=0;i<__hexes.size();++i) OUT << __hexes[i] << endl;
    // quads
    OUT << __quads.size() << endl;
    for (int i=0;i<__quads.size();++i) OUT << __quads[i] << endl;
    // lines
    OUT << __lines.size() << endl;
    for (int i=0;i<__lines.size();++i) OUT << __lines[i] << endl;

    // hanging & boundary
    this->write_data(__hanging_lines,OUT);
    this->write_data(__hanging_quads,OUT);
    this->write_data(__boundary_lines,OUT);
    this->write_data(__boundary_quads,OUT);
  
    OUT.close();

//    cerr << "write mesh:" << endl
//      << "\t" << __vertices.size() << " " << __hexes.size() << " " << __quads.size() << " " 
//      << __lines.size()<< endl;
//    cerr << "\t" << __hanging_lines.size() << " "
//      << __hanging_quads.size() << " "
//      << __boundary_lines.size() << " "
//      << __boundary_quads.size() << endl;
  }


  template<int DIM>
  void TriaContainer<DIM>::read_gip(const std::string& _fn) 
  {
    string fn = _fn;
    int name_size = fn.size();
    if(name_size<4) fn += ".gip";
    else if(fn.substr(name_size-4,4)!=".gip") fn += ".gip";

    ifstream IN(fn.c_str(),ios::in| ios::binary);
    if (!IN)
      {
	cerr << "TriaContainer::read_gip: could not open file " << fn << " for reading" << endl;
	abort();
      }
    reset();


    int tmp;
    char c;

    VERTEX tmpv;
    // vertices
    IN >> tmp; if (tmp>0) __vertices.resize(tmp);

    IN >> c;
    assert(c=='[');
    IN.read(reinterpret_cast<char*>(&__vertices[0][0]), sizeof(double) * DIM * __vertices.size());
    IN >> c;
    assert(c==']');

    // hexes
    IN >> tmp; if (tmp>0) __hexes.resize(tmp);
    for (int i=0;i<__hexes.size();++i) IN >> __hexes[i];
    // quads
    IN >> tmp; if (tmp>0) __quads.resize(tmp);
    for (int i=0;i<__quads.size();++i) IN >> __quads[i];
    // lines
    IN >> tmp; if (tmp>0) __lines.resize(tmp);
    for (int i=0;i<__lines.size();++i) IN >> __lines[i];
  

    // hanging & boundary
    this->read_data(__hanging_lines,IN);
    this->read_data(__hanging_quads,IN);
    this->read_data(__boundary_lines,IN);
    this->read_data(__boundary_quads,IN);
  
    assert(IN);
  
    IN.close();

//    cerr << "read mesh:" << endl
//      << "\t" << __vertices.size() << " " << __hexes.size() << " " << __quads.size() << " " 
//      << __lines.size()<< endl;
//    cerr << "\t" << __hanging_lines.size() << " "
//      << __hanging_quads.size() << " "
//      << __boundary_lines.size() << " "
//      << __boundary_quads.size() << endl;
  }






  template<>
  int TriaContainer<1>::first_smooth_flags(){ return 0; }

  template<>
  int TriaContainer<3>::first_smooth_flags(){ return 0; }

  //////////////////////////////////////////////////.  

  template<>
  int TriaContainer<2>::first_smooth_flags()
  {
    for (int l=0;l<nlines();++l)
      {
#define L line(l)
	if (L.type()!=0) continue;
	if (L.flag()==0) continue;
	if (line_at_boundary(L.id())) continue;
	if (L.slave()!=-1) continue;
	if (L.father()==-1) continue;

#define FL line(L.father())
	if (FL.flag()) continue;
#undef FL

#define Q quad(L.master())
	int lid;
	for (lid=0;lid<4;++lid) 
	  if (Q.line(lid)==l) break;
	assert(lid<4);
	// remove refinement flag
	if ((lid==0)||(lid==2))
	  Q.flag()^=REF_X;
	if ((lid==1)||(lid==3))
	  Q.flag()^=REF_Y;
	L.flag()=0;
#undef Q
	
#undef L
      }

  }
  









  template<>
  int TriaContainer<1>::smooth_flags(){ return 0; }

  template<>
  int TriaContainer<3>::smooth_flags(){ return 0; }

  //////////////////////////////////////////////////.  

  template<>
  int TriaContainer<2>::smooth_flags()
  {
    int smooth=0;
    

    // check, if left + right or bottom + top are flagged
    // but element itself is not in that direction

    for (int c=0;c<ncells();++c)
      {
#define Q quad(c)
	if ((line(Q.line(0)).flag())&&(line(Q.line(0)).flag())&&(!(Q.flag()&REF_X)))
	  {
	    if (set_refine_flag(c,REF_X))
	      {
		++smooth;
		assert(Q.flag()&REF_X);
	      }
	  }
	if ((line(Q.line(1)).flag())&&(line(Q.line(3)).flag())&&(!(Q.flag()&REF_Y)))
	  {
	    if (set_refine_flag(c,REF_Y))
	      {
		++smooth;
		assert(Q.flag()&REF_Y);
	      }
	  }
#undef Q
      }


    // do not allow hanging nodes which are to be resolved on the interface
    for (int c=0;c<ncells();++c)
      {
#define Q quad(c)
	if (Q.type()!=0) continue; // only active
	
	for (int l=0;l<4;++l)
	  {
#define L line(Q.line(l))
	    if (L.flag()==0) continue;
	    if (L.slave()==-1) continue;
	    
	    
	    int oi=0;
	    
	    for (int n=0;n<2;++n)
	      {
		const VERTEX& v = vertex(L.node(n));
		if ( ((fabs(v.y()-0.19)<1.e-8)||(fabs(v.y()-0.21)<1.e-8)) && 
		     (v.x()<=0.6+1.e-10)&&(v.x()>0.2) )
		  ++oi;
	      }
	    if (oi!=1) continue;

	    if (((l==0)||(l==2)) && (!(Q.flag()&REF_X)))
	      {
		++smooth;
		assert(set_refine_flag(c,REF_X));
	      }
	    if (((l==1)||(l==3)) && (!(Q.flag()&REF_Y)))
	      {
		++smooth;
		assert(set_refine_flag(c,REF_Y));
	      }
#undef L
	  }
	
#undef Q
	
      }
    




    return smooth;
  }
  


  //////////////////////////////////////////////////.

  template<>
  int TriaContainer<1>::resolve_flags(){ return 0; }
  
  //////////////////////////////////////////////////.  

  template<>
  int TriaContainer<2>::resolve_flags()
  {
    // check for hanging nodes on hanging nodes.
    // 
    // look for hanging lines, where child is flagged, but 
    // not the father.
    // resolve by flagging coarse quad
    int ac    = 0;
    int count = 0;

    do
      {
	count = 0;
	HASH_MAP::const_iterator hl;
	for (hl=__hanging_lines.begin();hl!=__hanging_lines.end();++hl)
	  {
#define L line(hl->second)
	    assert(L.type()==1);
	    if (L.flag()==1) continue;
	    if ((line(L.child(0)).flag()!=0)||(line(L.child(1)).flag()!=0))
	      {
		int lim = find_local_line_index_of_quad(L.master(), L.id());
		int lis = find_local_line_index_of_quad(L.slave(),  L.id());
		int rdm = ((lim==2)||(lim==0))?REF_X:REF_Y;
		int rds = ((lis==2)||(lis==0))?REF_X:REF_Y;
		
		if      (((quad(L.master()).type()&rdm)==0)&&
			 ((quad(L.master()).flag()&rdm)==0))
		  {
		    if (!set_refine_flag_on_quad(L.master()))
		      abort();
		  }
		else if (((quad(L.slave()). type()&rds)==0)&&
			 ((quad(L.slave()). flag()&rds)==0))
		  {
		    if (!set_refine_flag_on_quad(L.slave()))
		      abort();
		  }
		else abort();
		++count;
		++ac;
	      }
#undef L
	  }

	// check for hanging line on hanging line aniso
	if (!isotropic())
	  {
	    // list of new hanging lines...
	    set<int> new_hanging_lines;
	    for (int l=0;l<nlines();++l)
	      {
#define L line(l)
		if (L.type()!=0) continue;
		if (L.flag()==0) continue;
		if (L.master()==-1) continue;
		if (L.slave()==-1)  continue;
		
		int lim = find_local_line_index_of_quad(L.master(), l);
		int lis = find_local_line_index_of_quad(L.slave(),  l);
		int rdm = ((lim==2)||(lim==0))?REF_X:REF_Y;
		int rds = ((lis==2)||(lis==0))?REF_X:REF_Y;

		if (((quad(L.master()).flag()&rdm)==0)||
		    ((quad(L.slave()).flag() &rds)==0))
		  {
		    new_hanging_lines.insert(l);
		  }
		

#undef L
	      }
	    
	    HASH_MAP::const_iterator hlit;
	    for (hlit=__hanging_lines.begin();hlit!=__hanging_lines.end();++hlit)
	      new_hanging_lines.insert(hlit->second);
	    


	    set<int>::const_iterator nhl;
	    for (nhl=new_hanging_lines.begin();nhl!=new_hanging_lines.end();++nhl)
	      {
#define L line(*nhl)
		// not necessary... i think, but would fasten the program
		if (L.nchilds() != 0)
		  if (__hanging_lines.find(line(L.child(0)).node(0)) != __hanging_lines.end() 
		      || __hanging_lines.find(line(L.child(0)).node(1)) != __hanging_lines.end())
		    if (L.flag()==1) continue;
		//if (__hanging_lines.find(*nhl)!=__hanging_lines.end())
		//  if (L.flag()==1) continue;
		for (int n=0;n<2;++n)
		  {
		    if (__hanging_lines.find(L.node(n))!=__hanging_lines.end())
		      {
			HASH_MAP::const_iterator hll = __hanging_lines.find(L.node(n));

#define LL line(hll->second)
			int lim = find_local_line_index_of_quad(LL.master(), LL.id());
			int lis = find_local_line_index_of_quad(LL.slave(),  LL.id());
			int rdm = ((lim==2)||(lim==0))?REF_X:REF_Y;
			int rds = ((lis==2)||(lis==0))?REF_X:REF_Y;
			
			if      (((quad(LL.master()).type()&rdm)==0)&&
				 ((quad(LL.master()).flag()&rdm)==0))
			  {
			    if (!set_refine_flag_on_quad(LL.master()))
			      abort();
			    ++count;
			    ++ac;
			  }
			else if (((quad(LL.slave()). type()&rds)==0)&&
				 ((quad(LL.slave()). flag()&rds)==0))
			  {
			    if (!set_refine_flag_on_quad(LL.slave()))
			      abort();
			    ++count;
			    ++ac;
			  }
#undef LL
		      }
		  }
		

#undef L
	      }
	    
	  }
	


	// check flags
	if (isotropic())
	  for (int c=0;c<ncells();++c)
	    {
	      int fl = cell(c).flag();
	      if ((fl!=0)&&(fl!=3))
		{
		  cerr << "Flag " << fl << " on element " << c << endl;
		  abort();
		}
	    }
      }
    while (count != 0);
    return ac;
  }

  //////////////////////////////////////////////////.

  template<>
  int TriaContainer<3>::resolve_flags()
  {
    // check for hanging on hanging
    // for iso: 
    //   check active hexes. flag if quad or line 
    //     is refined and flagged for refinement.
    assert(isotropic());

    int count = 0;
    int ac    = 0;
     
    do
      {
	count = 0;
	for (int h=0;h<nhexes();++h)
	  {
#define H hex(h)
	    if (H.type()!=0) continue;
	    if (H.flag()!=0) continue;
	    // check quads
	    int q=0;
	    for (q=0;q<6;++q)
	      {
#define Q quad(H.quad(q))
		if (Q.type()==0) continue;
		// check childs (iso)
		assert(Q.type()==3);
		int c=0;
		for (c=0;c<4;++c)
		  if (quad(Q.child(c)).flag()!=0)
		    {
		      if (!set_refine_flag_on_hex(h))
			abort();
		      ++count; ++ ac;
		      break;
		    }
		if (c<4) break;
#undef Q
	      }
	    if (q<6) continue;
	    // check lines
	    int l=0;
	    for (l=0;l<12;++l)
	      {
#define L line(line_of_hex(h,l))
		if (L.type()==0) continue;
		int c=0;
		for (c=0;c<2;++c)
		  if (line(L.child(c)).flag()!=0)
		    {
		      if (!set_refine_flag_on_hex(h))
			abort();
		      ++count; ++ ac;
		      break;		      
		    }
		if (c<2) break;
#undef L
	      }
#undef H    
	  }
	
      }
    while (count!=0);
    
    return ac;
  }

  template class TriaContainer<2>;
  template class TriaContainer<3>;
}