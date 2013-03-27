/*----------------------------   element.h     ---------------------------*/
/*      $Id: element.h,v 1.4 2008/11/04 15:47:17 meidner Exp $                 */
#ifndef __element_H
#define __element_H
/*----------------------------   element.h     ---------------------------*/

#include <vector>
#include <cassert>
#include "mesharray.h"


/**
 *
 * Base class for an element with EDIM
 * This is mainly the hierarchy information, no link
 * to neighbors, ...
 *
 **/

#define REF_X 1
#define REF_Y 2
#define REF_Z 4
#define REF_ISO (REF_X|REF_Y|REF_Z)


namespace Tsuchimikado
{
  // EDIM: Type of Element. 2: quad, 3: hex  
  template <int EDIM>
    class Element
    {
#define N_MAX_CHILDS (EDIM*EDIM-EDIM+2)
#define N_SUBDATA    (2*EDIM)
    protected:

      // refinement type of element: 0 not refined, REF_X/Y/Z or combination
      int __type;
      // flags for refinement: as __type
      int __flag;
 

      // index of this element
      int               __id;
      // every element has one father
      int               __father;
      // reference to children, (depending on the dimension max 4 or 8)
      mesharray<N_MAX_CHILDS,int> __childs;
      // the nodes of the element. This should'nt be stored
      // in all elements in further versions.
      // see ...
      mesharray<N_MAX_CHILDS,int> __nodes;

      // in 2d: lines have one or two adjacent quads
      // in 3d: quad have one or two adjacent hexes 
      // !!!! lines in 2d do not need master/slave
      int __master;
      int __slave;
      // hexes have 6 quads, quads have 4 lines, lines have 2 ???
      mesharray<N_SUBDATA,int> __subdata;
  
    public:
  
      Element();
      // constructor to create a quad with reference to father
      // and to 4 nodes but without any other information
      Element(int f,int n1,int n2,int n3,int n4);

      // constructor to initialize the element with
      // father, childs and nodes
      Element(const int father,
	      const mesharray<N_MAX_CHILDS,int>& childs,
	      const mesharray<N_MAX_CHILDS,int>& nodes);

      //
      bool is_isotropic() const;

  
      // Access
      const int  id()     const { return __id; }
      const int type()   const { return __type; }
      const int flag()   const { return __flag; }
      const int  slave()  const { return __slave; }
      const int  master() const { return __master; }
  
      int&  id()     { return __id; }
      int&  type()   { return __type; }
      int&  flag()   { return __flag; }
      int&  slave()  { return __slave; }
      int&  master() { return __master; }
  

      const int n_max_childs() const{ return N_MAX_CHILDS; }
      const int nnodes() const      { return N_MAX_CHILDS; }
      const int nchilds() const
	{
	  if (__type==0) return 0;
	  if (EDIM==1) assert(__type==1);
	  if (EDIM==2) assert(__type==3);
	  
	  return n_max_childs();
	}
      const int father() const      { return __father; }
      const int child(int i) const  { assert(i<nchilds()); return __childs[i]; }
      const int node(int i) const   { assert(i<nnodes());  return __nodes[i]; }
      
      int& father()       { return __father; }
      int& child(int i)   { assert(i<nchilds()); return __childs[i]; }
      int& node(int i)    { assert(i<nnodes());  return __nodes[i]; }



      // Functions for the construction on an element. Will be removed from this class
      
      void init(const int father,
		const mesharray<N_MAX_CHILDS,int>& childs,
		const mesharray<N_MAX_CHILDS,int>& nodes);

      void new_by_lines(const Element<1>& l1,const Element<1>& l2,const Element<1>& l3,const Element<1>& l4);
      
      

      // Dimension specific access
      const int quad(int i) const 
	{ assert(EDIM==3); assert(i<6); return __subdata[i]; }
      const int line(int i) const 
	{ assert(EDIM==2); assert(i<4); return __subdata[i]; }
      int& quad(int i)
	{ assert(EDIM==3); assert(i<6); return __subdata[i]; }
      int& line(int i)
	{ assert(EDIM==2); assert(i<4); return __subdata[i]; }

  
      /**
       * print the elemendata
       **/
      void print() const;

      template<int X>
	friend std::ostream& operator<<(std::ostream& s, const Element<X>& E) ;
      template<int X>
	friend std::istream& operator>>(std::istream& s, Element<X>& E);
      
            
      // comparison
      // rotation, orientation does not matter
      // this must be specified for hexes
      const bool operator==(const Element<EDIM>& E) const;
#undef N_MAX_CHILDS
    };

}


/*----------------------------   element.h     ---------------------------*/
/* end of #ifndef __element_H */
#endif
/*----------------------------   element.h     ---------------------------*/
