#ifndef  __Q2Lps3d_h
#define  __Q2Lps3d_h

#include  "q23d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Lps3d

////
////
/////////////////////////////////////////////

class Q2Lps3d : public Q23d
{
public:

//
////  Con(De)structor 
//

  Q2Lps3d() : Q23d() {}
  ~Q2Lps3d() {}

  std::string GetName() const {return "Q2Lps3d";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
