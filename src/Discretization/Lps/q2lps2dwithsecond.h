#ifndef __Q2Lps2dWithSecond_h
#define __Q2Lps2dWithSecond_h

#include  "q2lps2d.h"
#include  "integratorlpswithsecond.h"

namespace Gascoigne
{

/**********************************************************/

class Q2Lps2dWithSecond : public Q2Lps2d
{
  protected:

  public:
   
    std::string GetName() const {return "Q2Lps2dWithSecond";}
    
    void BasicInit(const ParamFile* paramfile);
};

/**********************************************************/

}
#endif
