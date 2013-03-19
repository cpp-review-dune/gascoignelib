/*----------------------------   element.h     ---------------------------*/
/*      $Id: element.h,v 1.4 2008/11/04 15:47:17 meidner Exp $                 */
#ifndef __element_H
#define __element_H
/*----------------------------   element.h     ---------------------------*/

#include <vector>
#include <cassert>
#include "mesharray.h"

/**
 * a lot of work has to be done:
 *
 * - the data structure for a quad in 2 and 3d should be different. We store too much
 *   stuff here!!!
 *
 * - the access to the information is a big piece of shit! some functions, getting
 *   a line of a quad in 2d are functions of the element, get a line of a hex is a function
 *   in the triacontainer
 *
 * - 
 **/



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


namespace Tsuchimikado
{
/*   template<int EDIM> class Element; */
  
/*   template<int EDIM> std::ostream& operator<<(std::ostream& s, const Element<EDIM>& E); */
/*   template<int EDIM> std::istream& operator>>(std::istream& s,       Element<EDIM>& E); */

  
  template <int EDIM>
    class Element
    {
#define N_MAX_CHILDS (EDIM*EDIM-EDIM+2)
#define N_SUBDATA    (2*EDIM)
    protected:

      // refinement type of element
      int __type;
      // flags for refinement
      int __flag;
 

      // index of this element
      int               __id;
      // every element has one father
      int               __father;
      // depending on the dimension, 
      mesharray<N_MAX_CHILDS,int> __childs;
      // the nodes of the element. This should'nt be stored
      // in all elements in further versions.
      // see ...
      mesharray<N_MAX_CHILDS,int> __nodes;

      // perhaps this is not the correct place,
      // but all elements acting as boundary have a master and slave
      // we should have another class, where we know the overall dimension.
      // this could save a lot (or some) memory... later
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
	  if (EDIM==1)
	    {
	      if (__type==REF_X) return 2;
	      return 0;
	    }
	  else if (EDIM==2)
	    {
	      if ((__type==REF_X)||(__type==REF_Y)) return 2;
	      else if (__type==(REF_X|REF_Y)) return 4;
	      return 0;
	    }
	  else if (EDIM==3)
	    {
	      if ((__type==REF_X)||(__type==REF_Y)||(__type==REF_Z)) return 2;
	      else if ((__type==(REF_X|REF_Y))||(__type==(REF_X|REF_Z))||(__type==(REF_Y|REF_Z))) return 4;
	      else if (__type==(REF_X|REF_Y|REF_Z)) return 8;
	      return 0;
	    }
	  assert(0);
	}
      const int father() const      { return __father; }
      const int child(int i) const  { assert(i<nchilds()); return __childs[i]; }
      const int node(int i) const   { assert(i<nnodes()); return __nodes[i]; }


      int& father()       { return __father; }
      int& child(int i)   { assert(i<nchilds()); return __childs[i]; }
      int& node(int i)    { assert(i<nnodes()); return __nodes[i]; }

  
      // Construction
      void init(const int father,
		const mesharray<N_MAX_CHILDS,int>& childs,
		const mesharray<N_MAX_CHILDS,int>& nodes);

      void new_by_lines(const Element<1>& l1,
			const Element<1>& l2,
			const Element<1>& l3,
			const Element<1>& l4);
      
      

      // Dimension specific
      const int quad(int i) const 
	{ assert(EDIM==3); assert(i<6); return __subdata[i]; }
      const int line(int i) const 
	{
	  if (EDIM==2)
	    {
	      assert(i<4); return __subdata[i];
	    }
	  else assert(0);
	}
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
