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

  void BasicInit(const Edge*, int, const fixarray<2*DIM-2,int>&);
  void AddNodes(const Gascoigne::LocalVector&);

  const fixarray<2*DIM-2,int>&  GetVertex() const { return _vertex; }
  const Gascoigne::LocalVector& GetValue()  const { return _u; }
  const Edge&                   GetEdge()   const { return *_edge; }
  int                           GetCount()  const { return _count; }
  fixarray<2*DIM-2,double>      GetNorm()   const;

  void ShowStatistics() const;
};

/**********************************************************/

#endif
