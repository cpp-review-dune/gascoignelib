#ifndef  __index_h
#define  __index_h

#include  <map>
#include  <set>

#include  "stlio.h"
#include  "gascoigne.h"

/*------------------------------------------*/

class Index
{
 protected:

  typedef  std::map<int,int>  IntMap;

  Gascoigne::IntVector  vl2g, el2g, hl2g, ql2g;
  IntMap     vg2l, eg2l, hg2l, qg2l;

 public:

 protected:

  ////////////////
  // Local To Global
  ////////////////

  const Gascoigne::IntVector&  Vertexl2g() const { return vl2g;}
  const Gascoigne::IntVector&  Edgel2g()   const { return el2g;}
        Gascoigne::IntVector&  Edgel2g()         { return el2g;}

  int Edgel2g  (int i)  const { return el2g[i];}

 public:

  ////////////////
  // Constructors
  ////////////////

  Index();
  Index(const Index& I);
  Index& operator=(const Index& I);
  
  int Vertexl2g(int i)  const { return vl2g[i];}

        Gascoigne::IntVector&  Quadl2g()         { return ql2g;}
        Gascoigne::IntVector&  Hexl2g()          { return hl2g;}
  const Gascoigne::IntVector&  Hexl2g()    const { return hl2g;}
  const Gascoigne::IntVector&  Quadl2g()   const { return ql2g;}
  int Quadl2g  (int i)  const { return ql2g[i];}
  int Hexl2g   (int i)  const { return hl2g[i];}
  const IntMap&  Quadg2l()   const { return qg2l;}
  const IntMap&  Hexg2l()    const { return hg2l;}
  
  ////////////////
  // Sizes
  ////////////////

  int nnodes   ()    const  { return VertexSize();}
  int VertexSize ()   const { return vl2g.size();}
  int VertexGSize()   const { return vg2l.size();}
  int EdgeSize ()     const { return el2g.size();}
  int EdgeGSize()     const { return eg2l.size();}
  int HexSize ()      const { return hl2g.size();}
  int HexGSize()      const { return hg2l.size();}
  int QuadSize ()     const { return ql2g.size();}
  int QuadGSize()     const { return qg2l.size();}

  ////////////////
  // Global To Local
  ////////////////

  const IntMap&  Vertexg2l() const { return vg2l;}
  const IntMap&  Edgeg2l()   const { return eg2l;}
        IntMap&  Edgeg2l()         { return eg2l;}

  int Quadg2l  (int i) const
    {
      std::map<int,int>::const_iterator ip = qg2l .find(i);
      if(ip==qg2l.end()) {
	std::cerr << "Index:: Quadg2l" << std::endl;
	std::cerr << "there is no " << " "<< i << std::endl;
	assert(0);}
      return ip->second; 
    }
  int Hexg2l   (int i) const
    {
      std::map<int,int>::const_iterator ip = hg2l .find(i);
      if(ip==hg2l.end()) {
	std::cerr << "Index:: Hexg2l" << std::endl;
	std::cerr << "there is no " << " "<< i << std::endl;
	assert(0);}
      return ip->second; 
    }
  int Vertexg2l(int i) const
    {
      std::map<int,int>::const_iterator ip = vg2l .find(i);
      if(ip==vg2l.end()) {
	std::cerr << "Index:: Vertexg2l" << std::endl;
	std::cerr << "there is no " << " "<< i << std::endl;
	assert(0);}
      return ip->second; 
    }
  int Edgeg2l  (int i) const
    {
      std::map<int,int>::const_iterator ip = eg2l .find(i);
      if(ip==eg2l.end()) {
	std::cerr << "Index:: Edgeg2l" << std::endl;
	std::cerr << "there is no " << " "<< i << std::endl;
	assert(0);}
      return ip->second; 
    }

  ////////////////////////
  // Global To Local Check
  ////////////////////////

  // gibt -2 zurueck, falls globaler vertex nicht in levelmesh

  int Vertexg2lCheck(int i) const
  {
    IntMap::const_iterator ip = vg2l .find(i);
    if(ip==vg2l .end()) return -2;
    return ip->second;
  }
  int Edgeg2lCheck  (int i) const
  {
    IntMap::const_iterator ip = eg2l .find(i);
    if(ip==eg2l .end()) return -2;
    return ip->second;
  }
  int Quadg2lCheck  (int i) const
  {
    IntMap::const_iterator ip = qg2l .find(i);
    if(ip==qg2l .end()) return -2;
    return ip->second;
  }
  int Hexg2lCheck   (int i) const
  {
    IntMap::const_iterator ip = hg2l .find(i);
    if(ip==hg2l .end()) return -2;
    return ip->second;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Index& I);

  void  InitNodes (const Gascoigne::IntSet& nodes);
  void  InitEdges (const Gascoigne::IntSet& edges);
  void  InitQuads ();
  void  InitHexs  ();
};

#endif
