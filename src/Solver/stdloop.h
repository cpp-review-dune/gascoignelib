#ifndef  __StdLoop_h
#define  __StdLoop_h

#include  "extrapolator.h"
#include  "adaptordata.h"
#include  "basicloop.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

class StdLoop : public virtual BasicLoop
{

protected:

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_write;


  
  int _nmin, _nmax, _coarse;
  double _p;
  int    _random_coarsening;
  std::string _estimator, _extrapolate, _refiner;
  DoubleVector _JErr;
  Extrapolator    Extra;

  // new vectors

  DoubleVector ComputeFunctionals(VectorInterface& f, VectorInterface& u);

  DoubleVector GetExactValues() const;

  virtual void EtaVisu(std::string name, int i, const DoubleVector& eta) const;
  virtual void EtaCellVisu(std::string name, int i, const GlobalVector& eta) const;
  virtual void AdaptMesh(const DoubleVector& eta);
  virtual void AdaptMesh(const DoubleVector& eta,std::string refine_or_coarsen_step);
  virtual DoubleVector Functionals(VectorInterface& u, VectorInterface& f);
  virtual double Estimator(DoubleVector& eta, VectorInterface& u, VectorInterface& f);

public:

  StdLoop();
  ~StdLoop();

  void BasicInit(const ParamFile* paramfile,
		 const ProblemContainer* PC,
		 const FunctionalContainer* FC=NULL);

  void run(const std::string& problemlabel);
  void ClockOutput() const;
};
}
/*-----------------------------------------*/

#endif
