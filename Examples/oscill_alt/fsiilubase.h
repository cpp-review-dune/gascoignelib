/*----------------------------   fsiilubase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsiilubase_H
#define __fsiilubase_H
/*----------------------------   fsiilubase.h     ---------------------------*/


#include "gascoigne.h"

namespace Gascoigne
{
  
  
  class FSIIluBase
  {
  public:
    virtual void solve_fluid(GlobalVector& x) const=0;
    virtual void solve_solid(GlobalVector& x) const=0;

    virtual void modify_ilu_F_NS() { abort(); }
    virtual void modify_ilu_S() { abort(); }
    
    virtual void ConstructStructure(const IntVector& permlfuid, const IntVector& permsolid, const MatrixInterface& A)
    {
      std::cerr << "FSIIluBase::CS not written!" << std::endl;
      abort();
    }
    
  };


}



/*----------------------------   fsiilubase.h     ---------------------------*/
/* end of #ifndef __fsiilubase_H */
#endif
/*----------------------------   fsiilubase.h     ---------------------------*/
