#include "edge.h"

/*---------------------------------------------------*/

Edge& Edge::operator=(const Edge& e)
{
  c1 = e.master();
  c2 = e.slave();
  l1 = e.LocalMasterIndex();
  l2 = e.LocalSlaveIndex ();
}

/*---------------------------------------------------*/

std::pair<int,int> Edge::EdgeNeighbour(int i) const
{
  //const Edge&  E = edge(quad(i).edge(e));
  int in = master();
  int il = LocalMasterIndex();
  if(in==i)    
    {
      in = slave();
      il = LocalSlaveIndex();
    }
  return std::make_pair(in,il);
}

/*---------------------------------------------------*/

void Edge::swapping(int newindex)
{
  master() = newindex;
  LocalMasterIndex() = LocalSlaveIndex();
  LocalSlaveIndex()  = -1;
  slave() = -1;
}

/*---------------------------------------------------*/

void Edge::setmaster(int newindex, int newlocal)
{
  master() = newindex;
  LocalMasterIndex() = newlocal;
  LocalSlaveIndex()  = -1;
  slave() = -1;
}

/*---------------------------------------------------*/

std::ostream& operator<<(std::ostream &s, const Edge& A)
{
  s << A.master()  << " ";
  s << A.LocalMasterIndex()  << " ";
  s << A.slave()   << " ";
  s << A.LocalSlaveIndex() << " "<< std::endl;
  
  return s;
}

/*---------------------------------------------------*/

std::istream& operator>>(std::istream &s, Edge& A)
{
  s >> A.master() >> A.LocalMasterIndex() >> A.slave() >> A.LocalSlaveIndex();
}

/*---------------------------------------------------*/
