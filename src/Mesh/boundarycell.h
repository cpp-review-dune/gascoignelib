#ifndef __boundarycell_h
#define __boundarycell_h

#include  "cell.h"

#define NEDGES 1

template<int N>
class BoundaryCell : public Cell<N,NEDGES>
{
 protected:

  int mat, eiq, oq;

 public:

  BoundaryCell(int l = 0, int f = -1) : Cell<N,NEDGES>(l,f), mat(0), eiq(0), oq(0) {}
  BoundaryCell(const BoundaryCell& c) : Cell<N,NEDGES>(c)  , mat(c.material()), 
    eiq(c.edge_in_quad()), oq(c.of_quad()) {}

  int  nnchild ()     const { return N; }
  int  material()     const { return mat; }
  int& material()           { return mat; }
  int  edge_in_quad() const { return eiq; }
  int& edge_in_quad()       { return eiq; }
  int  of_quad()      const { return oq; }
  int& of_quad()            { return oq; }

  BoundaryCell<N>& operator=(const BoundaryCell<N>& c)
    {
      Cell<N,NEDGES>::operator=(c);
      mat = c.material();
      eiq = c.edge_in_quad();
      oq  = c.of_quad();
      return *this;
    }
  friend std::ostream& operator<<(std::ostream &s, const BoundaryCell& A)
    {
      s << " : ";
      s << A.vertex()  << " " << A.level()   << " ";
      s << A.father()  << " " << A.of_quad() << " ";
      s << A.edge_in_quad() << " @ ";
      s << A.nchilds() << " " << A.childs();
      
      
      /*s << "s-l-f-q  " << A.sleep() << " " << A.level();
	s << " " << A.father() << " " << A.of_quad() << std::endl;
	s << "Childs " << A.childs();
	s << std::endl;*/
      return s;
    }
  friend std::istream& operator>>(std::istream &s, BoundaryCell& A)
    {
      std::string symbol;
      int n;
      s >> symbol;
      if (symbol!=":")
	{
	  std::cout << "ERROR in BoundaryCell::operator>>" << std::endl;
	  exit(1);
	}
      s >> A.vertex() >> A.level();
      s >> A.father() >> A.of_quad() >> A.edge_in_quad();
      s >> symbol >> n;
      if (symbol!="@")
	{
	  std::cout << "ERROR in BoundaryCell::operator>>" << std::endl;
	  exit(1);
	}
      A.childs().resize(n);
      s >> A.childs();

      return s;
    }
};

#undef NEDGES

#endif
