#ifndef  __StdLoop_h
#define  __StdLoop_h

#include  "extrapolator.h"
#include  "adaptordata.h"
#include  "basicloop.h"

//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

namespace Gascoigne
{
class StdLoop : public BasicLoop
{
private:

  std::vector<const Functional*>   _FV;

protected:

  const std::vector<const Functional*>& GetFunctionals() const { return _FV;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_write;

  std::string _estimator, _extrapolate, _refiner;
  DoubleVector _JErr;
  Extrapolator    Extra;

  GlobalVector _GlobalErr;

  // new vectors

  DoubleVector ComputeFunctionals(MultiLevelGhostVector& f, MultiLevelGhostVector& u, const std::vector<const Functional*>& J) const;

  DoubleVector GetExactValues() const;

  virtual void EtaVisu(std::string name, int i, const DoubleVector& eta) const;
  virtual void AdaptMesh(const DoubleVector& eta);
  virtual DoubleVector Functionals(MultiLevelGhostVector& u, MultiLevelGhostVector& f);
  virtual double Estimator(DoubleVector& eta, MultiLevelGhostVector& u, MultiLevelGhostVector& f);

public:

  StdLoop();
  ~StdLoop();

  void BasicInit(const ParamFile* paramfile);

  void AddFunctional(const Functional* fv) { _FV.push_back(fv);}

  void run(const ProblemDescriptorInterface* PD);
  void ClockOutput() const;
};
}
/*-----------------------------------------*/

#endif
