#ifndef __Q23dWithSecond_h
#define __Q23dWithSecond_h

#include  "q23d.h"
#include  "integratorwithsecond.h"

/**********************************************************/

class Q23dWithSecond : public Gascoigne::Q23d
{
  protected:

  public:
   
    std::string GetName() const {return "Q23dWithSecond";}
    
    void BasicInit(const Gascoigne::ParamFile* paramfile);
};

/**********************************************************/

#endif
