#ifndef  __multigridmeshinterface_h
#define  __multigridmeshinterface_h

#include  "meshinterface.h"

/*--------------------------------------*/

class MultiGridMeshInterface
{
 public:
  
  virtual ~MultiGridMeshInterface() {}
  virtual int nlevels () const =0;
  virtual const MeshInterface& operator()(int l) const=0;
};

#endif
