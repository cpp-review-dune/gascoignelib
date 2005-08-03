#ifndef __edge_h
#define __edge_h

#include <iostream>
#include <utility>

/***********************************************/

namespace Gascoigne
{
class Edge
{
 protected:

  int c1, c2, l1, l2;  

 public:

  Edge()                { c1=c2=l1=l2=-1;}
  Edge(const Edge& e)   { *this=e; }
  Edge(int c, int l)    { c1=c,l1=l; c2=l2=-1;}

  int  master()           const { return c1;}
  int& master()                 { return c1;}
  int  slave ()           const { return c2;}
  int& slave ()                 { return c2;}

  int  LocalMasterIndex() const { return l1;}
  int& LocalMasterIndex()       { return l1;}
  int  LocalSlaveIndex () const { return l2;}
  int& LocalSlaveIndex ()       { return l2;}

  Edge& operator=(const Edge& e);

  void swapping(int);
  void setmaster(int,int);
  std::pair<int,int> EdgeNeighbour(int i) const;

  void BinWrite(std::ostream &s) const;
  void BinRead(std::istream &s);

  friend std::ostream& operator<<(std::ostream &s, const Edge& A);
  friend std::istream& operator>>(std::istream &s, Edge& A);
};
}

/***********************************************/

#endif
