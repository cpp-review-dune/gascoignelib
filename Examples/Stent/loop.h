/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "stdloop.h"

namespace Gascoigne
{
  
  class Loop : public StdLoop
  {
    
  public:
    
    void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,const FunctionalContainer* FC)
    {
      StdLoop::BasicInit(paramfile,PC,FC);
    }
    
    void run(const std::string& problemlabel);
  };
  
}


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
