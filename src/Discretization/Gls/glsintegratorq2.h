#ifndef  __GlsIntegratorQ2_h
#define  __GlsIntegratorQ2_h

#include  "glsintegrator.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GlsIntegratorQ2

////
////
/////////////////////////////////////////////

template<int DIM>
class GlsIntegratorQ2 : public GlsIntegrator<DIM>
{
public:

//
////  Con(De)structor 
//

  GlsIntegratorQ2<DIM>();
  ~GlsIntegratorQ2<DIM>() {}
};
}

#endif
