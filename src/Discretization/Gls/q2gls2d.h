#ifndef  __Q2Gls2d_h
#define  __Q2Gls2d_h

#include  "q22d.h"
#include  "q22dwithsecond.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Gls2d

////
////
/////////////////////////////////////////////

class Q2Gls2d : public Q22dWithSecond
{
public:

//
////  Con(De)structor 
//
  Q2Gls2d() : Q22dWithSecond() {}
  ~Q2Gls2d() {}

  std::string GetName() const {return "Q2Gls2d";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
