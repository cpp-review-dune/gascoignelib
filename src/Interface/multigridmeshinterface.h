#ifndef  __multigridmeshinterface_h
#define  __multigridmeshinterface_h

#include  "meshinterface.h"

/*--------------------------------------*/

namespace Gascoigne
{
  class MultiGridMeshInterface
  {
    private:

    protected:

    public:
      MultiGridMeshInterface() {}
      virtual ~MultiGridMeshInterface() {}

      virtual int nlevels () const =0;
      virtual const MeshInterface& operator()(int l) const=0;
  };
}

#endif
