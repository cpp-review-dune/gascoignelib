#ifndef __Q2Lps3dWithSecond_h
#define __Q2Lps3dWithSecond_h

#include  "q2lps3d.h"
#include  "integratorlpswithsecond.h"

/**********************************************************/

class Q2Lps3dWithSecond : public Gascoigne::Q2Lps3d
{
  protected:

  public:
   
    std::string GetName() const {return "Q2Lps3dWithSecond";}
    
    void BasicInit(const Gascoigne::ParamFile* paramfile);
};

/**********************************************************/

#endif
