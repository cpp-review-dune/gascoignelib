#ifndef  __PressureFilter_h
#define  __PressureFilter_h

#include "nvector.h"
#include "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PressureFilter : public DoubleVector 
{
 protected:

  IntVector component;
  double       domainsize;

 public:

  PressureFilter() ;
  ~PressureFilter();

  void SetComponents(const IntVector& c) { component = c;}
  bool Active() const { return component.size()>0;}

  void ReInit(int n);

  void AddDomainPiece(double val) { domainsize += val;}

  DoubleVector IntegrateVector(const GlobalVector& u) const;
  void SubtractMean(GlobalVector& u) const;
  void SubtractMeanAlgebraic(GlobalVector& u) const;
};
}

/*-----------------------------------------*/

#endif
