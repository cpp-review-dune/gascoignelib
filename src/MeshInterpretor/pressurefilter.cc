#include "pressurefilter.h"

/*-----------------------------------------*/

namespace Gascoigne
{
PressureFilter::PressureFilter() : nvector<double>(), domainsize(0.),  nhanging(0) {}

/*-----------------------------------------*/

PressureFilter::~PressureFilter() {}

/*-----------------------------------------*/

void PressureFilter::ReInit(int n, int nhn)
{
  resize(n);
  zero();
  domainsize = 0.;
  nhanging = nhn;
}

/*-----------------------------------------*/

DoubleVector PressureFilter::IntegrateVector(const GlobalVector& u) const
{
  assert(size());
  assert(Active());
  
  DoubleVector dst(u.ncomp(),0.);
  
  for (int j=0; j<u.n(); j++)
    {
      for (int i=0; i<component.size(); i++)
	{
	  int   c = component[i];
	  dst[c] += u(j,c)* (*this)[j];
	}      
    }
  return dst;
}

/*-----------------------------------------*/

void PressureFilter::SubtractMean(GlobalVector& u) const
{
  assert(size());
  assert(domainsize>0.);
  DoubleVector mean = IntegrateVector(u);
  
  for (int i=0; i<component.size(); i++)
    {
      int   comp = component[i];
      double sub = mean[comp]/domainsize;
      u.CompAdd(comp,-sub);
    }  
}

/*-----------------------------------------*/

void PressureFilter::SubtractMeanAlgebraic(GlobalVector& u) const
{
  for (int i=0; i<component.size(); i++)
    {
      int comp = component[i];
      double d = 0.;
      for (int j=0; j<u.n(); j++)
	{
	  d += u(j,comp)* (*this)[j];
	}      
      //d /= u.n();
      d /= u.n() - nhanging;
      u.CompAdd(comp,-d);
    }
}
}

/*-----------------------------------------*/

