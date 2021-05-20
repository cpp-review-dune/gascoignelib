/*----------------------------   dgbase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgbase_H
#define __dgbase_H
/*----------------------------   dgbase.h     ---------------------------*/

#include "baseq12d.h"
#include "baseq22d.h"

namespace Gascoigne {
class BASEQ12D : public BaseQ12d
{
public:
  static constexpr size_t N = 4;
};

class BASEQ22D : public BaseQ22d
{
public:
  static constexpr size_t N = 9;
};
} // namespace Gascoigne

/*----------------------------   dgbase.h     ---------------------------*/
/* end of #ifndef __dgbase_H */
#endif
/*----------------------------   dgbase.h     ---------------------------*/
