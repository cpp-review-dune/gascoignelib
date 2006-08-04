#ifndef __PiQ2_h
#define __PiQ2_h

#include "gascoignemesh.h"

namespace Gascoigne
{

/**********************************************************/

  class PiQ2
  {
    protected:
      const GascoigneMesh   *_MP;
      nvector<DoubleVector>  _q2weight;

    public:
      PiQ2();
      ~PiQ2() { }

      void Init(const MeshInterface* MI);
      void vmult(GlobalVector& y, const GlobalVector& x) const;
  };

/**********************************************************/

}

#endif
