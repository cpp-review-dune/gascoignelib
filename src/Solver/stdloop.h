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

class StdLoop : public BasicLoop
{
private:

  std::vector<const Functional*>   _FV;

protected:

  const std::vector<const Functional*>& GetFunctionals() const { return _FV;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_write;

  std::string _estimator, _extrapolate, _refiner;
  nvector<double> _JErr;
  Extrapolator    Extra;

  CompVector<double> _GlobalErr;

  // new vectors

  nvector<double> ComputeFunctionals(MultiLevelGhostVector& f, MultiLevelGhostVector& u, const std::vector<const Functional*>& J) const;

  nvector<double> GetExactValues() const;

  virtual void EtaVisu(std::string name, int i, const nvector<double>& eta) const;
  virtual void AdaptMesh(const nvector<double>& eta);
  virtual nvector<double> Functionals(MultiLevelGhostVector& u, MultiLevelGhostVector& f);
  virtual double Estimator(nvector<double>& eta, MultiLevelGhostVector& u, MultiLevelGhostVector& f);

public:

  StdLoop();
  ~StdLoop();

  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void AddFunctional(const Functional* fv) { _FV.push_back(fv);}

  void run(const ProblemDescriptorInterface* PD);
  void ClockOutput() const;
};

/*-----------------------------------------*/

#endif
