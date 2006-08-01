#ifndef __baseq22dPatchWithSecond_h
#define __baseq22dPatchWithSecond_h

#include  "baseq2patch.h"

namespace Gascoigne
{
/**************************************************/

  class BaseQ22dPatchWithSecond : public BaseQ22dPatch
  {
    public:
      BaseQ22dPatchWithSecond() : BaseQ22dPatch()
      {
        second = true;
      }
  };
}

#endif
