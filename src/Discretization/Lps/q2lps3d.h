#ifndef  __Q2Lps3d_h
#define  __Q2Lps3d_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Lps3d

////
////
/////////////////////////////////////////////

#include  "q23d.h"

namespace Gascoigne
{

class Q2Lps3d : public Q23d
{
public:

//
////  Con(De)structor 
//

  Q2Lps3d() : Q23d() {}
  ~Q2Lps3d() {}

  std::string GetName() const {return "Q2Lps";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
