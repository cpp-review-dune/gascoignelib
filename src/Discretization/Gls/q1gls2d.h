#ifndef  __Q1Gls2d_h
#define  __Q1Gls2d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments NewQ1Simple

////
////
/////////////////////////////////////////////

#include  "q12d.h"

namespace Gascoigne
{
class Q1Gls2d : public Q12d
{
protected:

public:

  //
  ////  Con(De)structor 
  //
  
  Q1Gls2d() : Q12d() {}
  ~Q1Gls2d() {}
  
  std::string GetName() const {return "Q1Gls2d";}
  
  void BasicInit(const ParamFile* pf);
};
}

#endif
