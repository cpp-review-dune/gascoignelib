#include  "index.h"

/*---------------------------------------------------*/

Index::Index() {}

/*---------------------------------------------------*/

Index::Index(const Index& I) 
{
  *this = I;
}

/*---------------------------------------------------*/

Index& Index::operator=(const Index& I)
{
  vl2g  = I.Vertexl2g();
  vg2l  = I.Vertexg2l() ;
  el2g  = I.Edgel2g();
  eg2l  = I.Edgeg2l();
  hl2g  = I.Hexl2g();
  hg2l  = I.Hexg2l();
  ql2g  = I.Quadl2g();
  qg2l  = I.Quadg2l();
}

/*---------------------------------------------------*/

std::ostream& operator<<(std::ostream& os, const Index& I)
{
  os << "Vertex l2g " << I.VertexSize()<<std::endl;
  os << I.Vertexl2g();
  os << "Vertex g2l " << I.VertexSize()<<std::endl;
  os << I.Vertexg2l() << " ";

  os << "Edge l2g " << I.EdgeSize()<<std::endl;
  os << I.Edgel2g();
  os << "Edge g2l " << I.EdgeSize()<<std::endl;
  os << I.Edgeg2l() << " ";

  os << "Hex l2g " << I.HexSize()<<std::endl;
  os << I.Hexl2g();
  os << "Hex g2l " << I.HexSize()<<std::endl;
  os << I.Hexg2l() << " ";

  os << "Quad l2g " << I.QuadSize()<<std::endl;
  os << I.Quadl2g();
  os << "Quad g2l " << I.QuadSize()<<std::endl;
  os << I.Quadg2l() << " ";

  return os;
}

/*---------------------------------------------------*/

void Index::InitNodes(const IntSet& nodes)
{
  typedef IntSet::const_iterator  iterator;
  int n = nodes.size();
  vl2g.memory(n);
  int count=0;
  iterator p=nodes.begin();
  for(int i=0;i<n;i++)
    {
      vl2g[count] = *p++;
      count++;
    }
  vg2l.clear();
  for(int i=0;i<n;i++)
    {
      vg2l.insert(std::make_pair(vl2g[i],i));
    }
}

/*---------------------------------------------------*/

void Index::InitEdges(const IntSet& edges)
{
  typedef IntSet::const_iterator  iterator;
  int  n = edges.size();
  el2g.memory(n);
  int count=0;
  for(iterator p=edges.begin();p!=edges.end();p++)
    {
      int i = *p;
      el2g[count++] = i;
    }
}

/*---------------------------------------------------*/

void Index::InitQuads()
{
  int  n  = ql2g.size();
  qg2l.clear();
  for(int i=0;i<n;i++)
    {
      qg2l.insert(std::make_pair(ql2g[i],i));
    }
}

/*---------------------------------------------------*/

void Index::InitHexs()
{
  int  n  = hl2g.size();
  hg2l.clear();
  for(int i=0;i<n;i++)
    {
      hg2l.insert(std::make_pair(hl2g[i],i));
    }
}

/*---------------------------------------------------*/
