#ifndef  __StdTimeLoop_h
#define  __StdTimeLoop_h

#include  "stdloop.h"
#include  "timeinfo.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class StdTimeLoop : public virtual StdLoop
{
protected:

  TimeInfo    _timeinfo;
  virtual std::string SolveTimePrimal(VectorInterface& u, VectorInterface& f);

  void TimeInfoBroadcast();
  void InitSolution(VectorInterface& u);
  void InitSolution(VectorInterface& u, VectorInterface& f);

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const ParamFile* paramfile);

  void run(const ProblemDescriptorInterface* PD);
  void adaptive_run(const ProblemDescriptorInterface* PD);
};
}

#endif
