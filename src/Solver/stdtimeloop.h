#ifndef  __StdTimeLoop_h
#define  __StdTimeLoop_h

#include  "stdloop.h"
#include  "timeinfo.h"

/*-----------------------------------------*/


class StdTimeLoop : public StdLoop
{
protected:

  TimeInfo    info;
  std::string SolveTimePrimal(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f, string name="Results/solve");

  void L2Projection(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f);
  void TimeInfoBroadcast();

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const std::string& pfile, const ProblemDescriptorInterface* PD);

  void run();
  void adaptive_run();
};


#endif
