#ifndef  __StdTimeLoop_h
#define  __StdTimeLoop_h

#include  "stdloop.h"
#include  "timeinfo.h"

/*-----------------------------------------*/


class StdTimeLoop : public StdLoop
{
protected:

  TimeInfo    info;
  virtual std::string SolveTimePrimal(MultiLevelGhostVector& u, MultiLevelGhostVector& f, std::string name="Results/solve");

  void TimeInfoBroadcast();

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const Gascoigne::ParamFile* paramfile);
  void InitSolution(MultiLevelGhostVector& u, MultiLevelGhostVector& f);

  void run(const ProblemDescriptorInterface* PD);
  void adaptive_run(const ProblemDescriptorInterface* PD);
};


#endif
