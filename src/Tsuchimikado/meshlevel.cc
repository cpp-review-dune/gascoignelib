#include "meshlevel.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace Tsuchimikado
{

  
  //////////////////////////////////////////////////
  
  template<>
  int MeshLevel<2>::size() const { return nquads(); }
  template<>
  int MeshLevel<3>::size() const { return nhexes(); }

  template<>
    MeshLevel<2>::CIT MeshLevel<2>::begin() const { return __quads.begin(); }
  template<>
    MeshLevel<3>::CIT MeshLevel<3>::begin() const { return __hexes.begin(); }
  template<>
    MeshLevel<2>::IT MeshLevel<2>::begin() { return __quads.begin(); }
  template<>
    MeshLevel<3>::IT MeshLevel<3>::begin() { return __hexes.begin(); }

  template<>
    MeshLevel<2>::CIT MeshLevel<2>::end() const { return __quads.end(); }
  template<>
    MeshLevel<3>::CIT MeshLevel<3>::end() const { return __hexes.end(); }
  template<>
    MeshLevel<2>::IT MeshLevel<2>::end() { return __quads.end(); }
  template<>
    MeshLevel<3>::IT MeshLevel<3>::end() { return __hexes.end(); }
      

  template<>
    const int& MeshLevel<2>::operator[](size_t n) const { return __quads[n]; }
  template<>
    const int& MeshLevel<3>::operator[](size_t n) const { return __hexes[n]; }
  template<>
    int& MeshLevel<2>::operator[](size_t n) { return __quads[n]; }
  template<>
    int& MeshLevel<3>::operator[](size_t n) { return __hexes[n]; }

  template<>
    void MeshLevel<2>::push_back(const int& e) {__quads.push_back(e);  }
  template<>
    void MeshLevel<3>::push_back(const int& e) {__hexes.push_back(e);  }
  
  



  // ==================================================

  template<int DIM>
  void MeshLevel<DIM>::init_active()
  {
    assert(DIM==__TC->dimension());
  
    clear();
    
    for (int c=0;c<__TC->ncells();++c)
      if (__TC->cell(c).nchilds()==0) this->push_back(c);
    
    for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;

    post_init();
    init_hanging();
    if (!mesh_ok())
      abort();
  }

  // ==================================================

  template<int DIM>
  void MeshLevel<DIM>::init_from_meshlevel(const MeshLevel<DIM>& ML, const HASH_SET& dont_coarse)
  {
    // set of elements which are coarsened
    __changed.clear();

    clear();
    // temporary set for new cells
    set<int> newcells;

    // hash set of cells in fine mesh for fast access
    HASH_SET oldcells;
    for (CIT it = ML.begin();it!=ML.end();++it)
      oldcells.insert(*it);

    // add new cells
    for (CIT it = ML.begin();it!=ML.end();++it)
      {
	// if all siblings are in fine mesh, add father, 
	// otherwise add cell itself
	int ci = *it;
	const CELL& C = __TC->cell(ci);
	int fa = C.father();
	if ((fa!=-1) && (dont_coarse.find(ci)==dont_coarse.end()) )
	  {
	    const CELL& F = __TC->cell(fa);
	    int index=0;
	    for (index=0;index<F.nchilds();++index)
	      {
		int ch = F.child(index);
		if (oldcells.find(ch)==oldcells.end()) break;
	      }
	    // add father or cell itself
	    if (index==F.nchilds())
	      {
		newcells.insert(fa);
		__changed.insert(fa);
	      }
	    else 
	      newcells.insert(ci);
	  }
	else newcells.insert(ci);
      }
    // convert to vector
    for (set<int>::const_iterator i=newcells.begin();i!=newcells.end();++i)
      push_back(*i);

    __Cg2l.clear();
    for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;



    init_hanging();
    
    bool ok;
    do
      {
	ok = repair_mesh();
	init_hanging();
      }
    while (!ok);

    post_init();
  }

  template<int DIM>
  void MeshLevel<DIM>::post_init()
  {
    ////////////// for 3d, init quads
    set<int> tmp;

    if (DIM==3)
      {
	__quads.clear(); tmp.clear();
	for (CIT it=begin();it!=end();++it)
	  {
	    int cell = this->operator[](*it);
	    for (int nn=0;nn<6;++nn)
	      tmp.insert(__TC->cell(cell).quad(nn));
	  }
	std::copy(tmp.begin(), tmp.end(), std::back_inserter(__quads));
      }
    
    ///////////// init lines (2d and 3d)
    __lines.clear(); tmp.clear();
    for (int i=0;i<__quads.size();++i)
      {
	int quad = __quads[i];
	for (int nn=0;nn<4;++nn)
	  tmp.insert(__TC->quad(quad).line(nn));
      }
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(__lines));    
    

    ///////////////// init nodes (2d and 3d)
    __nodes.clear(); tmp.clear();
    for (CIT it=begin();it!=end();++it)
      {
	int cell = *it;
	for (int nn=0;nn<__TC->cell(cell).nnodes();++nn)
	  tmp.insert(__TC->cell(cell).node(nn));
      }
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(__nodes));

    // std::cerr << "mesh: nodes/lines/quads/hexes:  "
    // 	      << __nodes.size() << " " 
    // 	      << __lines.size() << " " 
    // 	      << __quads.size() << " " 
    // 	      << __hexes.size() << std::endl;
  }
  

  // ==================================================

  template<>
  void MeshLevel<2>::init_hanging()
  {
    __hanging_lines.clear();
    __hanging_quads.clear();


    // a line is a hanging line, if 
    //   it belongs to a quad in the mesh
    //   it is not at the boundary
    //   a child of the neighboring quad is in the same mesh
    for (int c=0;c<this->size();++c)
      {
	int cell = Celll2g(c);
	for (int l=0;l<4;++l)
	  {
	    int line = __TC->cell(cell).line(l);
	    if (__TC->line_at_boundary(line)) continue; // at boundary?
	    int othercell = __TC->line(line).master();
	    if (othercell==cell) othercell=__TC->line(line).slave();
	    if (othercell==-1) continue;                // line is fine
	    
	    if (__TC->cell(othercell).type()==0) continue; // othercell not refined
	    assert(__TC->cell(othercell).nchilds()==4);
	    for (int oc=0;oc<4;++oc)
	      if (CellInMesh(__TC->cell(othercell).child(oc)))
		__hanging_lines.insert(line);
	  }
      }

    __hanging_nodes_on_lines.clear();
    __hanging_nodes_on_quads.clear();
    for (HASH_SET::const_iterator it=__hanging_lines.begin();
	 it!=__hanging_lines.end();++it)
      __hanging_nodes_on_lines[__TC->middle_node(__TC->line(*it))]=*it;
    
  }

  // ==================================================

  template<>
  void MeshLevel<3>::init_hanging()
  {
    assert(__TC->isotropic());
    
    __hanging_lines.clear();
    __hanging_quads.clear();
    // create set of all lines and quads in this mesh
    HASH_SET lines,quads;
    for (CIT it=begin();it!=end();++it)
      {
	for (int q=0;q<6;++q)
	  quads.insert(__TC->cell(*it).quad(q));
	for (int l=0;l<12;++l)
	  lines.insert(__TC->line_of_hex(*it,l));
      }
    
    // lines/quads hang, if childs are in same mesh.
    HASH_SET::const_iterator it;
    for (it=lines.begin();it!=lines.end();++it)
      {
	int li = *it;
	if (__TC->line(li).nchilds()==0) continue;
	int c0 = __TC->line(li).child(0);
	int c1 = __TC->line(li).child(1);
	if (lines.find(c0)!=lines.end()) __hanging_lines.insert(li);
	if (lines.find(c1)!=lines.end()) __hanging_lines.insert(li);
      }
    for (it=quads.begin();it!=quads.end();++it)
      {
	int li = *it;
	if (__TC->quad(li).nchilds()==0) continue;

	assert(__TC->quad(li).type()==3);

	for (int ch=0;ch<4;++ch)
	  if (quads.find(__TC->quad(li).child(ch))!=quads.end()) __hanging_quads.insert(li);
      }

    __hanging_nodes_on_lines.clear();
    __hanging_nodes_on_quads.clear();
    for (HASH_SET::const_iterator it=__hanging_lines.begin();
	 it!=__hanging_lines.end();++it)
      __hanging_nodes_on_lines[__TC->middle_node(__TC->line(*it))]=*it;
    for (HASH_SET::const_iterator it=__hanging_quads.begin();
	 it!=__hanging_quads.end();++it)
      __hanging_nodes_on_quads[__TC->middle_node(__TC->quad(*it))]=*it; 
  }

  // --------------------------------------------------

  template<>
  bool MeshLevel<2>::mesh_ok()
  {
    // every line must:
    //   - appear twice
    //   - once and on the boundary
    //   - once and in hanging
    //   - once and father in hanging
    HASH_MAP lc;
    for (int c=0;c<size();++c)
      for (int l=0;l<4;++l)
	lc[__TC->cell(Celll2g(c)).line(l)]++;

    for (HASH_MAP::const_iterator it = lc.begin();it!=lc.end();++it)
      {
	if (it->second==2) continue;
	int line = it->first;
	if (__TC->line_at_boundary(line)) continue;
	if (line_is_hanging(line)) continue;
	
	int father = __TC->line(line).father();
	assert(father!=-1);
	if (line_is_hanging(father)) continue;
	return false;
      }
    return true;
  }

  // ----------------------------------------

  template<>
  bool MeshLevel<3>::mesh_ok()
  {
    abort();
//     // mesh has a problem, if hanging line/quad hangs on hanging line/quad.
//     __double_hanging_nodes.clear();
//     HASH_SET::const_iterator it;
//     for (it=__hanging_lines.begin();it!=__hanging_lines.end();++it)
//       {
// 	for (int n=0;n<2;++n)
// 	  {
// 	    int node = __TC->line(*it).node(n);
// 	    if (__hanging_nodes_on_lines.find(node)!=__hanging_nodes_on_lines.end())
// 	      __double_hanging_nodes.insert(node);
// 	    if (__hanging_nodes_on_quads.find(node)!=__hanging_nodes_on_quads.end())
// 	      __double_hanging_nodes.insert(node);
// 	  }
//       }
//     for (it=__hanging_quads.begin();it!=__hanging_quads.end();++it)
//       {
// 	for (int n=0;n<4;++n)
// 	  {
// 	    int node = __TC->quad(*it).node(n);
// 	    if (__hanging_nodes_on_lines.find(node)!=__hanging_nodes_on_lines.end())
// 	      __double_hanging_nodes.insert(node);
// 	    if (__hanging_nodes_on_quads.find(node)!=__hanging_nodes_on_quads.end())
// 	      __double_hanging_nodes.insert(node);
// 	  }
//       }
//     return (__double_hanging_nodes.size()==0);
  }

  // --------------------------------------------------

  template<>
  bool MeshLevel<2>::repair_mesh()
  {
    // there may be no double-hanging lines, e.g.
    //    lines in the mesh with father-father also in the mesh
    HASH_MAP lc;
    for (int c=0;c<size();++c)
      for (int l=0;l<4;++l)
	lc[__TC->cell(Celll2g(c)).line(l)]++;

    set<int> takebackcoarsening;

    for (HASH_MAP::const_iterator it = lc.begin();it!=lc.end();++it)
      {
	int line = it->first;

	if (it->second==2) continue;                // line twice     => ok
	if (__TC->line_at_boundary(line)) continue; // at boundary    => ok
	if (line_is_hanging(line)) continue;        // hanging        => ok
	int father = __TC->line(line).father();
	if (line_is_hanging(father)) continue;      // father hanging => ok

	
	// PROBLEM! 2 possibilities:
	//    * father-father of line is in the mesh => do nothing wait for other case
	//    * child-child of line is in the mesh   => fix master/slave of line
	
	
	// check if father-father is in the mesh
	int fatherfather = -1;
	if (father!=-1) fatherfather = __TC->line(father).father();
	if (lc.find(fatherfather)!=lc.end()) continue; // ok, wait for the other case

	// Check, that at least one child-child is in the mesh
	bool ok = false;
	assert(__TC->line(line).type()==1);
	for (int c=0;c<2;++c)
	  {
	    int child = __TC->line(line).child(c);
	    assert(child!=-1);
	    if (__TC->line(child).type()!=0)
	      for (int cc=0;cc<2;++cc)
		if (lc.find(__TC->line(child).child(cc))!=lc.end()) ok=true;
	  }
	if (!ok) 
	  {
	    cerr << "Mesh is broken!";
	    print_gnuplot2("meshes");
	    abort();
	  }
	
	int dontcoarse = __TC->line(line).master();
	assert(dontcoarse!=-1);
	if (!CellInMesh(dontcoarse)) dontcoarse = __TC->line(line).slave();
	assert(dontcoarse!=-1);
	assert(CellInMesh(dontcoarse));
	
	takebackcoarsening.insert(dontcoarse);
      }
     
    if (takebackcoarsening.size()==0) return true;
    for (set<int>::const_iterator it = takebackcoarsening.begin(); 
	 it!=takebackcoarsening.end();++it)
      {
	int cell = *it;
	assert(CellInMesh(cell));
	assert(__TC->cell(cell).type()==3);
	
	int oldnr = Cellg2l(cell);  // replace old cell by the four childs
	(*this)[oldnr]=__TC->cell(cell).child(0);
	for (int c=1;c<4;++c)
	  push_back(__TC->cell(cell).child(c));
      }

    __Cg2l.clear();
    for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;
    
    return false;


    
//     // double hanging lines are not ok.
//     for (HASH_SET::const_iterator it=__hanging_lines.begin();
// 	 it!=__hanging_lines.end();++it)
//       {
// 	int line = *it;
// 	if (__TC->line(line).type()==0) continue; // line is fine
// 	for (int l=0;l<2;++l)
// 	  if (line_is_hanging(__TC->line(line).child(0)))
// 	    {
// 	      abort();
	      
// 	    }
//       }
	   
//     // resolve double hanging nodes by taking back coarsening.
//     assert(__double_hanging_nodes.size()>0);

//     HASH_SET problem_lines;
    
    
//     for (HASH_MAP::const_iterator it=__double_hanging_nodes.begin();it!=__double_hanging_nodes.end();++it)
//       {
// 	//cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX " << it->first << endl;
// 	assert(__hanging_nodes_on_lines.find(it->first)!=__hanging_nodes_on_lines.end());
// 	problem_lines.insert(__hanging_nodes_on_lines[it->first]);
// 	problem_lines.insert(__hanging_nodes_on_lines[it->second]);
//       }    // we need to find the elements which have the problem_lines as boundary lines.
//     // these elements must be replaced by finer ones.
//     // they are changed elements
 
//     // um zu ueberpruefenn, wie viel elemente gerade in problem_lines enthalten sind
//     //cout << "wir haben jetzt " << problem_lines.size() << " doppel haengende knoten, naemlich " 
//     //	 <<  *problem_lines.begin() << endl;   

// //     if (problem_lines.size() >= 1)
// //       {
// // 	cout << "Vertex bei "<< *problem_lines.begin() << " mit " 
// // 	     << problem_lines.size() << endl;
// // 	//for (it=problem_lines.begin(); it!=problem_lines.end(); it++)
// // 	//  cout << __TC->vertex(*it) << endl;
// //       }
// //     cout << endl;

//     set<int> add;
//     //in changed are all elements, which shall be coarsened for coarser grid
//     for (HASH_SET::const_iterator it = __changed.begin();it!=__changed.end();++it)
//       for (int l=0;l<4;++l)
// 	{
// 	  if (problem_lines.find(__TC->quad(*it).line(l))!=problem_lines.end())
// 	    {
// 	      add.insert(*it);
// 	    }
// 	}
//     assert(add.size()>0);
//     for (set<int>::const_iterator i1 = add.begin();i1!=add.end();++i1)
//       {
// 	//	assert(__TC->quad(*i1).type()==3);
// 	__changed.erase(*i1);
// 	assert(__Cg2l.find(*i1)!=__Cg2l.end());
// 	int li = __Cg2l[*i1];
// 	this->operator[](li) = __TC->quad(*i1).child(0);
// 	for (int c=1;c<__TC->quad(*i1).nchilds();++c)
// 	  push_back(__TC->quad(*i1).child(c));
//       }
//     __Cg2l.clear();
//     for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;
  }

  // --------------------------------------------------

  template<>
  bool MeshLevel<3>::repair_mesh()
  {
    abort();
//     assert(__double_hanging_nodes.size()>0);
//     HASH_SET problem_lines,problem_quads;

//     HASH_SET::const_iterator it;
//     for (it=__double_hanging_nodes.begin();it!=__double_hanging_nodes.end();++it)
//       {
// 	if (__hanging_nodes_on_lines.find(*it)!=__hanging_nodes_on_lines.end())
// 	  problem_lines.insert(__hanging_nodes_on_lines[*it]);
// 	else if (__hanging_nodes_on_quads.find(*it)!=__hanging_nodes_on_quads.end())
// 	  problem_quads.insert(__hanging_nodes_on_quads[*it]);
// 	else assert(0);
//       }

//     // we need to find the elements which have the problem_lines as boundary lines.
//     // these elements must be replaced by finer ones.
//     // they are changed elements
//     set<int> add;
//     for (it = __changed.begin();it!=__changed.end();++it)
//       {
// 	for (int l=0;l<12;++l)
// 	  if (problem_lines.find(__TC->line_of_hex(*it,l))!=problem_lines.end())
// 	    add.insert(*it);
// 	for (int q=0;q<6;++q)
// 	  if (problem_quads.find(__TC->hex(*it).quad(q))!=problem_quads.end())
// 	    add.insert(*it);
//       }
//     assert(add.size()>0);

//     for (set<int>::const_iterator i1 = add.begin();i1!=add.end();++i1)
//       {
// 	assert(__TC->hex(*i1).type()==7);
// 	__changed.erase(*i1);
// 	assert(__Cg2l.find(*i1)!=__Cg2l.end());
// 	int li = __Cg2l[*i1];
// 	this->operator[](li) = __TC->hex(*i1).child(0);
// 	for (int c=1;c<8;++c)
// 	  push_back(__TC->hex(*i1).child(c));
//       }
//     __Cg2l.clear();
//     for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;
  }

  template<int DIM>
  void MeshLevel<DIM>::print_gnuplot(const std::string& fname) const
  {
    assert(DIM==2);
    ofstream OUT(fname.c_str());
    for (CIT it=begin();it!=end();++it)
      {
	for (int l=0;l<4;++l)
	  {
	    Gascoigne::Vertex<DIM> v1 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(0));
	    Gascoigne::Vertex<DIM> v2 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(1));
	    OUT << v1 << endl << v2 << endl << endl;
	  }
      }
    
    OUT.close();
  }

  template<int DIM>
  void MeshLevel<DIM>::print_gnuplot2(std::string fname) const
  {
    assert(DIM==2);

    stringstream data;
    data << fname << "_data";

    ofstream OUT(data.str().c_str());
    for (CIT it=begin();it!=end();++it)
      {
	for (int l=0;l<4;++l)
	  {
	    Gascoigne::Vertex<DIM> v1 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(0));
	    Gascoigne::Vertex<DIM> v2 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(1));
	    OUT << v1 << endl << v2 << endl << endl;
	  }
      }
    OUT.close();

    double minx=0,maxx=0,miny=0,maxy=0;
    minx=maxx=__TC->vertex(0).x();
    miny=maxy=__TC->vertex(0).y();


    for (int i=0;i<__TC->nvertices();++i)
      {
	minx = std::min(minx, __TC->vertex(i).x());
	maxx = std::max(maxx, __TC->vertex(i).x());
	miny = std::min(miny, __TC->vertex(i).y());
	maxy = std::max(maxy, __TC->vertex(i).y());
      }

    double dx = maxx-minx;    assert(dx>0);
    double dy = maxx-minx;    assert(dy>0);

    OUT.open(fname.c_str());
    OUT << "set xrange [" << minx-0.05*dx << ":" << maxx+0.05*dx << "]" << endl;
    OUT << "set yrange [" << miny-0.05*dy << ":" << maxy+0.05*dy << "]" << endl;
    OUT << "set sty da li" << endl;

    OUT << "plot \"" << data.str() << "\"" << endl;
    OUT << "pause -1" << endl;
    OUT.close();
  }




  // ==================================================
  template<int DIM>
    MeshLevel<DIM>::MeshLevel() : __TC(0)
    {
      abort();
    }

  template<int DIM>
    MeshLevel<DIM>::MeshLevel(const TriaContainer<DIM>& TC) : __TC(&TC)
    {}
  

  template<int DIM>
    void MeshLevel<DIM>::clear() 
    {
      __nodes.clear();
      __lines.clear();
      __quads.clear();
      __hexes.clear();
      __Cg2l.clear();
      __hanging_lines.clear();
      __hanging_quads.clear();
      __hanging_nodes_on_lines.clear();
      __hanging_nodes_on_quads.clear();
      __double_hanging_nodes.clear();
      __changed.clear();
    }

  template class MeshLevel<2>;
  template class MeshLevel<3>;
 
}

