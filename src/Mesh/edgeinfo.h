#ifndef __EdgeInfo_h
#define __EdgeInfo_h

#include "edge.h"
#include "gascoigne.h"

/**********************************************************/

template<int DIM>
class EdgeInfo
{
 protected:

  int                    _count;
  fixarray<2*DIM-2,int>  _vertex;
  Gascoigne::LocalVector _u;
  const Edge*            _edge;

 public:

  EdgeInfo<DIM>() {}
  ~EdgeInfo<DIM>() {}

  void basicInit(const Edge*, int, const fixarray<2*DIM-2,int>&);
  void addNodes(const Gascoigne::LocalVector&);
  void showStatistics() const;

  const fixarray<2*DIM-2,int>& vertex() const;

  int getCount() const;
  const Gascoigne::LocalVector& getValue() const;
  fixarray<2*DIM-2,double> getNorm() const;
  const Edge& getEdge() const;
};

/**********************************************************/

template EdgeInfo<2>;
template EdgeInfo<3>;

#endif
