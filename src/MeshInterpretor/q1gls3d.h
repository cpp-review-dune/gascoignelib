#ifndef  __Q1Gls3d_h
#define  __Q1Gls3d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Gls3d

////
////
/////////////////////////////////////////////

#include  "q13d.h"
#include  "glsintegrator.h"

namespace Gascoigne
{
class Q1Gls3d : public Q13d
{
protected:

/*   GlsIntegrator<3> GlsInt;   */

public:

  //
  ////  Con(De)structor 
  //
  
  Q1Gls3d() : Q13d() {}
  ~Q1Gls3d() {}
  
  std::string GetName() const {return "Q1Gls3d";}
  
  void BasicInit(const ParamFile* pf);
};
}

#endif
