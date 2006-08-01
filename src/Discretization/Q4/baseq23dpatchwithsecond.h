#ifndef __baseq23dPatchWithSecond_h
#define __baseq23dPatchWithSecond_h

#include  "baseq23dpatch.h"

namespace Gascoigne
{
/**************************************************/

  class BaseQ23dPatchWithSecond : public BaseQ23dPatch
  {
    public:
      BaseQ23dPatchWithSecond() : BaseQ23dPatch()
      {
        second = true;
      }
  };
}

#endif
