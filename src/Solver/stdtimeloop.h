#ifndef  __StdTimeLoop_h
#define  __StdTimeLoop_h

#include  "stdloop.h"
#include  "timeinfo.h"

/*-----------------------------------------*/


class StdTimeLoop : public StdLoop
{
protected:

  TimeInfo    _timeinfo;
  virtual std::string SolveTimePrimal(MultiLevelGhostVector& u, MultiLevelGhostVector& f);

  void TimeInfoBroadcast();
  void InitSolution(MultiLevelGhostVector& u);
  void InitSolution(MultiLevelGhostVector& u, MultiLevelGhostVector& f);

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void run(const ProblemDescriptorInterface* PD);
  void adaptive_run(const ProblemDescriptorInterface* PD);
};


#endif
