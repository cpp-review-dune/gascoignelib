/*----------------------------   cudaloop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef cudaloop_H
#define cudaloop_H
/*----------------------------   cudaloop.h     ---------------------------*/

#include <string>

#include "stdloop.h"

namespace Gascoigne {
class FunctionalContainer;
class ParamFile;
class ProblemContainer;

/**
 * Loop to add very simple time stepping
 **/
class CudaLoop : public StdLoop
{
public:
  void BasicInit(const ParamFile&,
                 const ProblemContainer*,
                 const FunctionalContainer*) override;

  // void timerun(const std::string&);
};

} // namespace Gascoigne

/*----------------------------   cudaloop.h     ---------------------------*/
/* end of #ifndef cudaloop_H */
#endif
/*----------------------------   cudaloop.h     ---------------------------*/
