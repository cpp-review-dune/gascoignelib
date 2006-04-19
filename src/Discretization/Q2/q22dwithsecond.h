#ifndef __Q22dWithSecond_h
#define __Q22dWithSecond_h

#include  "q22d.h"
#include  "integratorwithsecond.h"

namespace Gascoigne
{

/**********************************************************/

class Q22dWithSecond : public virtual Q22d
{
  protected:

  public:
   
    std::string GetName() const {return "Q22dWithSecond";}
    
    void BasicInit(const ParamFile* paramfile);
};

/**********************************************************/

}
#endif
