#ifndef  __LocalLoop_h
#define  __LocalLoop_h

#include  "stdloop.h"

#include  "meshagent.h"

/*-----------------------------------------*/

namespace Gascoigne
{

  class LocalLoop : public StdLoop
    {
    public:
      void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC=NULL)
      {
        GetMeshAgentPointer() = new MeshAgent;

        GetMeshAgent()->AddPeriodicMapping(90,92, new StdPeriodicMapping);
        GetMeshAgent()->AddPeriodicMapping(92,90, new StdPeriodicMapping);

        StdLoop::BasicInit(paramfile, PC, FC);
      };

    };
  
}

/*-----------------------------------------*/

#endif // __LocalLoop_h
