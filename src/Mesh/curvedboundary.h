#ifndef  __CurvedBoundary_h
#define  __CurvedBoundary_h

#include  "boundaryfunction.h"
#include  "vertex.h"

/*---------------------------------------------*/

class SinusBoundary :  public virtual BoundaryFunction<2>
{
public :
  double operator()(const Vertex2d& c) const;
};

#endif
