#ifndef __Q2Lps3dWithSecond_h
#define __Q2Lps3dWithSecond_h

#include  "q2lps3d.h"
#include  "integratorlpswithsecond.h"

namespace Gascoigne
{

/**********************************************************/

class Q2Lps3dWithSecond : public Q2Lps3d
{
  protected:

  public:
   
    std::string GetName() const {return "Q2Lps3dWithSecond";}
    
    void BasicInit(const ParamFile* paramfile);
};

/**********************************************************/

}
#endif
