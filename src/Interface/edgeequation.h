#ifndef __edgeequation_h
#define __edgeequation_h

//
#include "gascoigne.h"
#include "vertex.h"
#include "feminterface.h"

namespace Gascoigne
{
  class EdgeEquation
  {
  public:
    EdgeEquation(){}
    virtual ~EdgeEquation(){}

    // initialize an integration point on the edge. U1 and U2 are the solutions from the two sides, n is the
    // normal vector as seen from element 1.
    virtual void point_edge(double h, const FemFunction& U1, const FemFunction& U2, const Vertex2d& v, const Vertex2d& n) const{}
    virtual void point_edge(double h, const FemFunction& U1, const FemFunction& U2, const Vertex3d& v, const Vertex3d& n) const{}

    // integration on one edge
    // test-functions from both sides are given at the same time, to be sorted into
    // b1 and b2.
    virtual void EdgeForm(VectorIterator b1, VectorIterator b2,
			  const FemFunction& U1, const FemFunction& U2,
			  const TestFunction& N1, const TestFunction& N2) const
    { std::cerr << "EdgeEquation::EdgeForm not written" << std::endl; abort(); }

    
  };
  
}

#endif
