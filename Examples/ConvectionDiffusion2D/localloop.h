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
#include  "problemdescriptor1.h"

class LocalLoop : public StdLoop
{
protected:

  ProblemDescriptor1 LPD;

public:


//
////  Con(De)structor 
//
  LocalLoop() : StdLoop() {}
  ~LocalLoop() {}

  void BasicInit(const ParamFile* paramfile);
  void run();
};


#endif
