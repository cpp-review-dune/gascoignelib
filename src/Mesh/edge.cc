#include "edge.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
Edge& Edge::operator=(const Edge& e)
{
  c1 = e.master();
  c2 = e.slave();
  l1 = e.LocalMasterIndex();
  l2 = e.LocalSlaveIndex ();

  return *this;
}

/*---------------------------------------------------*/

pair<int,int> Edge::EdgeNeighbour(int i) const
{
  //const Edge&  E = edge(quad(i).edge(e));
  int in = master();
  int il = LocalMasterIndex();
  if(in==i)    
    {
      in = slave();
      il = LocalSlaveIndex();
    }
  return make_pair(in,il);
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

ostream& operator<<(ostream &s, const Edge& A)
{
  s << A.master()  << " ";
  s << A.LocalMasterIndex()  << " ";
  s << A.slave()   << " ";
  s << A.LocalSlaveIndex() << " "<< endl;
  
  return s;
}

/*---------------------------------------------------*/

istream& operator>>(istream &s, Edge& A)
{
  s >> A.master() >> A.LocalMasterIndex() >> A.slave() >> A.LocalSlaveIndex();

  return s;
}
}

/*---------------------------------------------------*/
