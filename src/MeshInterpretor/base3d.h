#ifndef __base3d_h
#define __base3d_h

#include  "base.h"

/**************************************************/

class Base3d : public Base
{
 protected:

  mutable nvector<double>             N;
  mutable std::vector<Vertex3d>       DN;

  mutable Vertex3d  bn, bt;

 public:
  
  Base3d()  {}
  const Vertex3d&  normal () const {return bn;}
  const Vertex3d&  tangent() const {return bt;}
};

#endif
