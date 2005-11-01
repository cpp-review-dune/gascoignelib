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

  virtual void TimeInfoBroadcast();
  void InitSolution(VectorInterface& u);
  void InitSolution(VectorInterface& u, VectorInterface& f);

public:

  StdTimeLoop() : StdLoop() {}

  void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC,
		 const FunctionalContainer* FC=NULL);

  void run(const std::string& problemlabel);
  void adaptive_run(const std::string& problemlabel);
};
}

#endif
