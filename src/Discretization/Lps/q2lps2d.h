#ifndef  __Q2Lps2d_h
#define  __Q2Lps2d_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Lps2d

////
////
/////////////////////////////////////////////

#include  "q22d.h"

namespace Gascoigne
{

class Q2Lps2d : public Q22d
{

public:

//
////  Con(De)structor 
//

  Q2Lps2d() : Q22d() {}
  ~Q2Lps2d() {}

  std::string GetName() const {return "Q2Lps2d";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
