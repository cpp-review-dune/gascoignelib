#ifndef  __PressureFilter_h
#define  __PressureFilter_h

#include "nvector.h"
#include "gascoigne.h"

/*-----------------------------------------*/

class PressureFilter : public nvector<double> 
{
 protected:

  nvector<int> component;
  double       domainsize;

 public:

  PressureFilter() ;
  ~PressureFilter();

  void SetComponents(const nvector<int>& c) { component = c;}
  bool Active() const { return component.size()>0;}

  void ReInit(int n);

  void AddDomainPiece(double val) { domainsize += val;}

  nvector<double> IntegrateVector(const Gascoigne::GlobalVector& u) const;
  void SubtractMean(Gascoigne::GlobalVector& u) const;
  void SubtractMeanAlgebraic(Gascoigne::GlobalVector& u) const;
};

/*-----------------------------------------*/

#endif
