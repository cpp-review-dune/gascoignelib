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

  nvector<double> ComputeFunctionals(NewMultiLevelGhostVector& f, NewMultiLevelGhostVector& u, const std::vector<const Functional*>& J) const;

  nvector<double> GetExactValues() const;

  virtual void EtaVisu(std::string name, int i, const nvector<double>& eta);
  virtual void AdaptMesh(const nvector<double>& eta);
  virtual nvector<double> Functionals(NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f);
  virtual double Estimator(nvector<double>& eta, NewMultiLevelGhostVector& u, NewMultiLevelGhostVector& f);

public:

  StdLoop();
  ~StdLoop();

  void BasicInit(const Gascoigne::ParamFile* paramfile);

  void SetFunctionals(const std::vector<const Functional*>& fv) { _FV =  fv;}

  void run(const ProblemDescriptorInterface* PD);
  void ClockOutput() const;
};

/*-----------------------------------------*/

#endif
