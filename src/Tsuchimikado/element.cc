#include "element.h"

using namespace std;



namespace Tsuchimikado
{

  // Constructors

  //  -----------------------------

  template<int DIM>
  Element<DIM>::Element() :
    __type(0), __flag(0), __id(-1),
    __father(-1), __master(-1), __slave(-1)		
  {
    __childs   = -1;
    __nodes    = -1;
    __subdata  = -1;
  }

  // ------------------------------

  template<int DIM>
  Element<DIM>::Element(int f,int n1,int n2,int n3,int n4):
    __type(0), __flag(0), __id(-1),
    __father(-1), __master(-1), __slave(-1)		
  {
    assert(DIM==2);
    __nodes[0]=n1;
    __nodes[1]=n2;
    __nodes[2]=n3;
    __nodes[3]=n4;
  }

  // **************************************************

  template<int DIM>
  Element<DIM>::Element(const int father,
			const mesharray<DIM*DIM-DIM+2,int>& childs,
			const mesharray<DIM*DIM-DIM+2,int>& nodes) :
    __type(0), __flag(0), __id(-1),
    __father(-1), __master(-1), __slave(-1)		
  { 
    init (father, childs, nodes); 
  }
  
  // **************************************************

  template<int EDIM>
  bool Element<EDIM>::is_isotropic() const
  {
    if (nchilds()==0)              return true;
    if (EDIM==1) return (__type==1);
    if (EDIM==2) return (__type==3);
    if (EDIM==3) return (__type==7);
    abort();
  }
  
  // **************************************************

  template<int EDIM>
  void Element<EDIM>::init(const int father,
			   const mesharray<EDIM*EDIM-EDIM+2,int>& childs,
			   const mesharray<EDIM*EDIM-EDIM+2,int>& nodes)
  {
    __father   = father;
    __childs   = childs;
    __nodes    = nodes;
  }

  // **************************************************

  template<int EDIM>
  void Element<EDIM>::print() const
  {
    if (__id==0)
      std::cout << "EDIM" << "\t"
		<< "id" << "\t"
		<< "type" << "\t"
		<< "flags" << "\t"
		<< "master" << "\t"
		<< "slave" << "\t"
		<< "father" << "\t"
		<< "nchilds" << "\t"
		<< "childs" << "\t"
		<< "subdata" << "\t"
		<< "nodes" << endl;
    
    std::cout << EDIM << "\t"
	      << __id << "\t"
	      << __type << "\t"
	      << __flag << "\t"
	      << __master << "\t"
	      << __slave << "\t"
	      << __father << "\t"
	      << nchilds() << "\t"
	      << __childs << "\t"
	      << __subdata << "\t"
	      << __nodes << endl;
  }


  template<>
  const bool Element<1>::operator==(const Element<1>& E) const
  {
    if ((node(0)==E.node(0))&&(node(1)==E.node(1))) return true;
    if ((node(1)==E.node(0))&&(node(0)==E.node(1))) return true;
    return false;
  }

  template<>
  const bool Element<2>::operator==(const Element<2>& E) const
  {
    // that's not good!

    // we have to check, if in any rotation and in
    // any of the two directions, the nodes fall together

    // check 4 rotations, forward
    for (int rot=0;rot<4;++rot)
      {
	int i;
	for (i=0;i<4;++i) if (E.node(i)!=node((i+rot)%4)) break;
	if (i==4) return true;
      }
    // check 4 rotations, backward
    for (int rot=0;rot<4;++rot)
      {
	int i;
	for (i=0;i<4;++i) if (E.node(3-i)!=node((i+rot)%4)) break;
	if (i==4) return true;
      }
    return false;
  }

  template<>
  const bool Element<3>::operator==(const Element<3>& E) const
  { // no need for this function
    abort();
  }

  template<int EDIM>
  std::ostream& operator<<(std::ostream& s, const Element<EDIM>& E) 
  {
    s << E.__type << " "
      << E.__flag << " "
      << E.__id << " "
      << E.__father << " "
      << E.__childs << " "
      << E.__nodes << " "
      << E.__master << " "
      << E.__slave << " "
      << E.__subdata << endl;
    return s;
  }
  template<int EDIM>
  std::istream& operator>>(std::istream& s, Element<EDIM>& E)
  {
    s >> E.__type
      >> E.__flag
      >> E.__id
      >> E.__father
      >> E.__childs
      >> E.__nodes
      >> E.__master
      >> E.__slave 
      >> E.__subdata;
    return s;
  }



  // **************************************************


  template<int DIM>
  void Element<DIM>::new_by_lines(const Element<1>& l1,
				  const Element<1>& l2,
				  const Element<1>& l3,
				  const Element<1>& l4)
  {
    assert(DIM==2);

    this->__father   = -1;
    this->__childs   = -1;
    this->__nodes    = -1;

    if ((l1.node(1)==l2.node(0))||(l1.node(1)==l2.node(1)))
      {
	this->__nodes[0]=l1.node(0);
	this->__nodes[1]=l1.node(1);
      }
    else if ((l1.node(0)==l2.node(0))||(l1.node(0)==l2.node(1)))
      {
	this->__nodes[0]=l1.node(1);
	this->__nodes[1]=l1.node(0);
      }
    else assert(0);

    if (this->__nodes[1]==l2.node(0))      this->__nodes[2]=l2.node(1);
    else if (this->__nodes[1]==l2.node(1)) this->__nodes[2]=l2.node(0);
    else assert(0);

    if (this->__nodes[2]==l3.node(0))      this->__nodes[3]=l3.node(1);
    else if (this->__nodes[2]==l3.node(1)) this->__nodes[3]=l3.node(0);
    else assert(0);
    assert(((l4.node(0)==this->__nodes[3])&&(l4.node(1)==this->__nodes[0]))||
	   ((l4.node(1)==this->__nodes[3])&&(l4.node(0)==this->__nodes[0])));
    
    for (int i=0;i<DIM*DIM-DIM+2;++i) assert(this->__nodes[i]!=-1);    

    this->line(0)=l1.id();
    this->line(1)=l2.id();
    this->line(2)=l3.id();
    this->line(3)=l4.id();
  }


  // --------------------------------------------------

  template std::ostream& operator<<(std::ostream& s, const Element<1>& E);
  template std::ostream& operator<<(std::ostream& s, const Element<2>& E);
  template std::ostream& operator<<(std::ostream& s, const Element<3>& E);
  template std::istream& operator>>(std::istream& s, Element<1>& E);
  template std::istream& operator>>(std::istream& s, Element<2>& E);
  template std::istream& operator>>(std::istream& s, Element<3>& E);



  template class Element<1>;
  template class Element<2>;
  template class Element<3>;

}
