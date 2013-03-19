#include "meshlevel.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace Tsuchimikado
{
  

  // ==================================================

  template<int DIM>
  void MeshLevel<DIM>::init_active()
  {
    assert(DIM==__TC->dimension());
  
    clear();
    __changed.clear();
    

    for (int c=0;c<__TC->ncells();++c)
      if (__TC->cell(c).nchilds()==0) this->push_back(c);
    
    __Cg2l.clear();
    for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;

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
    for (IT it = ML.begin();it!=ML.end();++it)
      oldcells.insert(*it);

    // add new cells
    for (IT it = ML.begin();it!=ML.end();++it)
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
    while (!mesh_ok())
      {
	repair_mesh();
	init_hanging();
      }
  }

  // ==================================================

  template<>
  void MeshLevel<2>::init_hanging()
  {
    __hanging_lines.clear();
    __hanging_quads.clear();

    // count appearance of all lines, not at the boundary.
    map<int,int> lc;
    for (IT it=begin();it!=end();++it)
      for (int l=0;l<4;++l)
	{
	  int line = __TC->cell(*it).line(l);
	  if (__TC->line_at_boundary(line)) continue;
	  if(lc.find(line)==lc.end())
	    lc[line] = 1;
	  else lc[line]++;
	}
    
    // line is a regular hanging line if
    // - it appears only once
    // - and one child is in mesh
    //
    // line is a double hanging line if
    // - it appears only once
    // - children are not in mesh
    // - some grandchildren are
    
    for (map<int,int>::const_iterator it = lc.begin();it!=lc.end();++it)
      {
	if (it->second==2) continue;
	assert(it->second==1);
	int line = it->first;
	if (__TC->line(line).type()==0) continue;
	int c[2] = {__TC->line(line).child(0),__TC->line(line).child(1)};	
	
	for (int j=0;j<2;++j)
	  {
	    // regular hanging
	    if (lc.find(c[j])!=lc.end()) __hanging_lines.insert(line);
	    // check for double hanging
	    if (__TC->line(c[j]).type()!=0)
	      {
		int gc[2] = {__TC->line(c[j]).child(0),__TC->line(c[j]).child(1)};
		for (int k=0;k<2;++k)
		  if (lc.find(gc[k])!=lc.end()) 
		    {
		      __hanging_lines.insert(line);
		      __hanging_lines.insert(c[j]);
		    }
	      }
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
    for (IT it=begin();it!=end();++it)
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
    // mesh has a problem, if hanging line hangs on hanging line.
    __double_hanging_nodes.clear();
    HASH_SET::const_iterator it;
    for (it=__hanging_lines.begin();it!=__hanging_lines.end();++it)
      {
	for (int n=0;n<2;++n)
	  {
	    int node = __TC->line(*it).node(n);
	    if (__hanging_nodes_on_lines.find(node)!=__hanging_nodes_on_lines.end())
	      {
		int middle_node = __TC->middle_node(__TC->line(*it));
		if (__double_hanging_nodes.find(node)!=__double_hanging_nodes.end())
		  {
		    // ???



// 		    print_gnuplot2("Z");
// 		    cout << __TC->vertex(node) << endl;
// 		    cerr << " double_hanging_node " << node << " 2mal!" << endl;
// 		    abort();
		  }
		__double_hanging_nodes[node]=middle_node;
	      }
	  }
      }

//     if (__double_hanging_nodes.size() > 0)
//       cout<< "soviel hab ich jetzt: " << __double_hanging_nodes.size() 
// 	  << " und das ist der  Knoten: "<< __double_hanging_nodes.begin()->first << endl 
// 	  << " und das ist der 2te  Knoten: "<< __double_hanging_nodes.begin()->second << endl 
// 	  << " der erste hat den vertex: " << __TC->vertex(__double_hanging_nodes.begin()->first) << endl
// 	  << " der zweite hat den vertex: " << __TC->vertex(__double_hanging_nodes.begin()->second) << endl;

    //Ausgabe des Zwischengitters
    print_gnuplot2("Y");
 

    return (__double_hanging_nodes.size()==0);
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
  void MeshLevel<2>::repair_mesh()
  {
    // resolve double hanging nodes by taking back coarsening.
    assert(__double_hanging_nodes.size()>0);

    HASH_SET problem_lines;
    
    
    for (HASH_MAP::const_iterator it=__double_hanging_nodes.begin();it!=__double_hanging_nodes.end();++it)
      {
	//cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXX " << it->first << endl;
	assert(__hanging_nodes_on_lines.find(it->first)!=__hanging_nodes_on_lines.end());
	problem_lines.insert(__hanging_nodes_on_lines[it->first]);
	problem_lines.insert(__hanging_nodes_on_lines[it->second]);
      }    // we need to find the elements which have the problem_lines as boundary lines.
    // these elements must be replaced by finer ones.
    // they are changed elements
 
    // um zu ueberpruefenn, wie viel elemente gerade in problem_lines enthalten sind
    //cout << "wir haben jetzt " << problem_lines.size() << " doppel haengende knoten, naemlich " 
    //	 <<  *problem_lines.begin() << endl;   

//     if (problem_lines.size() >= 1)
//       {
// 	cout << "Vertex bei "<< *problem_lines.begin() << " mit " 
// 	     << problem_lines.size() << endl;
// 	//for (it=problem_lines.begin(); it!=problem_lines.end(); it++)
// 	//  cout << __TC->vertex(*it) << endl;
//       }
//     cout << endl;

    set<int> add;
    //in changed are all elements, which shall be coarsened for coarser grid
    for (HASH_SET::const_iterator it = __changed.begin();it!=__changed.end();++it)
      for (int l=0;l<4;++l)
	{
	  if (problem_lines.find(__TC->quad(*it).line(l))!=problem_lines.end())
	    {
	      add.insert(*it);
	    }
	}
    assert(add.size()>0);
    for (set<int>::const_iterator i1 = add.begin();i1!=add.end();++i1)
      {
	//	assert(__TC->quad(*i1).type()==3);
	__changed.erase(*i1);
	assert(__Cg2l.find(*i1)!=__Cg2l.end());
	int li = __Cg2l[*i1];
	this->operator[](li) = __TC->quad(*i1).child(0);
	for (int c=1;c<__TC->quad(*i1).nchilds();++c)
	  push_back(__TC->quad(*i1).child(c));
      }
    __Cg2l.clear();
    for (int c=0;c<this->size();++c) __Cg2l[(*this)[c]]=c;
  }

  // --------------------------------------------------

  template<>
  void MeshLevel<3>::repair_mesh()
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
    for (IT it=begin();it!=end();++it)
      {
	for (int l=0;l<4;++l)
	  {
	    MeshVertex<DIM> v1 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(0));
	    MeshVertex<DIM> v2 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(1));
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
    for (IT it=begin();it!=end();++it)
      {
	for (int l=0;l<4;++l)
	  {
	    MeshVertex<DIM> v1 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(0));
	    MeshVertex<DIM> v2 = __TC->vertex(__TC->line(__TC->quad(*it).line(l)).node(1));
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

  template class MeshLevel<2>;
  template class MeshLevel<3>;
 
}

