#ifndef __Q2Lps2dWithSecond_h
#define __Q2Lps2dWithSecond_h

#include  "q2lps2d.h"
#include  "integratorlpswithsecond.h"

/**********************************************************/

class Q2Lps2dWithSecond : public Gascoigne::Q2Lps2d
{
  protected:

  public:
   
    std::string GetName() const {return "Q2Lps2dWithSecond";}
    
    void BasicInit(const Gascoigne::ParamFile* paramfile);
};

/**********************************************************/

#endif
