#ifndef __cell_h
#define __cell_h

#include  "fixarray.h"
#include  <string>
#include  "gascoigne.h"

/*-----------------------------------------------------------*/

// N = no of nodes
// E = no of edges/faces

/*-----------------------------------------------------------*/

namespace Gascoigne
{
template<int N, int E>
class Cell :  public fixarray<N,int>   /* das sind die vertex-no. */
{
 protected:

  /* Data */

  int              qlevel, qfather;
  IntVector         qchilds;
  fixarray<E,int>  qedges;            /* edge numbers */

 public:

  /* Constructors */

  Cell() : 
    fixarray<N,int>(), 
    qlevel(0), 
    qfather(-1) 
    { 
      qedges=-1; 
    }

  Cell(const Cell& c) : 
    fixarray<N,int>(c) , 
    qlevel(c.level()), 
    qfather(c.father()), 
    qchilds(c.childs()), 
    qedges(c.edges()) {}

  Cell(int l, int f) : 
    fixarray<N,int>(-17), 
    qlevel(l), 
    qfather(f) {}

  /* Operators */

  Cell<N,E>& operator=(const Cell<N,E>& c)
    {
      vertex() = c.vertex();
      qlevel   = c.level(); 
      qfather  = c.father(); 
      qchilds.memory(c.nchilds());
      qchilds  = c.childs();
      qedges   = c.edges();
      return *this;
    }
  bool operator==(const Cell<N,E>& c) const
    {
      if (vertex()==c.vertex()) return 1;
      return 0;
    }

  /* Zugriff */

  int   level  ()     const { return qlevel;}          
  int&  level  ()           { return qlevel;}          
  int   father ()     const { return qfather;}          
  int&  father ()           { return qfather;}          
  bool  sleep  ()     const { return qchilds.size()!=0;}          
  int   nchilds()     const { return qchilds.size();}          
  int   nvertexs()    const { return N;}          
  int   child (int i) const { return qchilds[i];}          
  int&  child (int i)       { return qchilds[i];}          
  int   vertex(int i) const { return (*this)[i];}          
  int&  vertex(int i)       { return (*this)[i];}          
  int   edge  (int i) const { return qedges[i];}          
  int&  edge  (int i)       { return qedges[i];}          

  const IntVector&  childs()       const { return qchilds;}          
        IntVector&  childs()             { return qchilds;}          
  const fixarray<N,int>& vertex() const { return (*this);}
        fixarray<N,int>& vertex()       { return (*this);}
  const fixarray<E,int>& edges()  const { return qedges;}
        fixarray<E,int>& edges()        { return qedges;}

  /* Functions */

  template<int M>
  void vertex_loc2glob(fixarray<M,int>& ig, 
		 const fixarray<M,int>& il) const
    {
      typename fixarray<M,int>::iterator        gp=ig.begin();
      typename fixarray<M,int>::const_iterator  lp=il.begin();
      while(lp!=il.end())  *gp++ = (*this)[*lp++];
    }

  int global2local(int gi) const;
  
  friend std::ostream& operator<<(std::ostream &s, const Cell& A)
    {
      s << A.vertex()  << " ";
      s << A.level()   << " ";
      s << A.father()  << " @ ";
      s << A.nchilds() << " " << A.childs();
      s << " : " << A.edges();
      s << std::endl;
      
      return s;
    }
  
  friend std::istream& operator>>(std::istream& s, Cell& A) 
    {
      std::string symbol;
      int  n;
      s >> A.vertex();
      s >> A.level() ;
      s >> A.father();
      s >> symbol;
      if (symbol!="@")
	{
	  std::cout << "ERROR: Cell::operator>>" << std::endl;
	  exit(1);
	}
      s >> n;
      A.childs().resize(n);
      s >> A.childs();
      s >> symbol;
      if (symbol!=":")
	{
	  std::cout << "ERROR: Cell::operator>>" << std::endl;
	  exit(1);
	}
      s >> A.edges();
      return s;
    }
};

template<int N, int E>
inline int Cell<N,E>::global2local(int gi) const
{
  for (int i=0; i<N; i++)
    {
      if (vertex(i)==gi) return i;
    }
  return -1;
}
}

/*---------------------------------------------------*/

#endif
