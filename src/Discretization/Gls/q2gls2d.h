#ifndef  __Q2Gls2d_h
#define  __Q2Gls2d_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Gls2d

////
////
/////////////////////////////////////////////

#include  "q22d.h"

namespace Gascoigne
{

class Q2Gls2d : public Q22d
{
public:

//
////  Con(De)structor 
//
  Q2Gls2d() : Q22d() {}
  ~Q2Gls2d() {}

  std::string GetName() const {return "Q2Gls";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
