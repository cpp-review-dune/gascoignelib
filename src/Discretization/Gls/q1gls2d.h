#ifndef  __Q1Gls2d_h
#define  __Q1Gls2d_h

#include  "q12d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments NewQ1Simple

////
////
/////////////////////////////////////////////

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
