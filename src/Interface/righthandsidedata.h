#ifndef  __righthandsidedata_h
#define  __righthandsidedata_h

#include  "gascoigne.h"
#include  "vertex.h"
#include  <string>
#include  "nvector.h"
#include  "application.h"

//////////////////////////////////////////////
///
///@brief
/// Interface class for RightHandSideData

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

namespace Gascoigne
{
class RightHandSideData : public Application
{
protected:

public:

  RightHandSideData() : Application() {}
  ~RightHandSideData() {}
};
}

#endif
