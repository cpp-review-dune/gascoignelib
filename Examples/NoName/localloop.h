#ifndef  __LocalLoop_h
#define  __LocalLoop_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalLoop

////
////
/////////////////////////////////////////////

#include  "stdloop.h"


class LocalLoop : public StdLoop
{
private:


protected:


public:


//
////  Con(De)structor 
//

  LocalLoop() : StdLoop() {}
  ~LocalLoop() {}

  void BasicInit(const Gascoigne::ParamFile* pfile);
};


#endif
