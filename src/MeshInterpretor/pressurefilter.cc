#include "pressurefilter.h"

/*-----------------------------------------*/

PressureFilter::PressureFilter() : nvector<double>(), domainsize(0.) {}

/*-----------------------------------------*/

PressureFilter::~PressureFilter() {}

/*-----------------------------------------*/

void PressureFilter::ReInit(int n)
{
  resize(n);
  zero();
  domainsize = 0.;
}

/*-----------------------------------------*/

nvector<double> PressureFilter::IntegrateVector(const Gascoigne::GlobalVector& u) const
{
  assert(size());
  assert(Active());
  
  nvector<double> dst(u.ncomp(),0.);
  
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

void PressureFilter::SubtractMean(Gascoigne::GlobalVector& u) const
{
  assert(size());
  assert(domainsize>0.);
  nvector<double> mean = IntegrateVector(u);
  
  for (int i=0; i<component.size(); i++)
    {
      int   comp = component[i];
      double sub = mean[comp]/domainsize;
      u.CompAdd(comp,-sub);
    }  
}

/*-----------------------------------------*/

void PressureFilter::SubtractMeanAlgebraic(Gascoigne::GlobalVector& u) const
{
  for (int i=0; i<component.size(); i++)
    {
      int comp = component[i];
      double d = 0.;
      for (int j=0; j<u.n(); j++)
	{
	  d += u(j,comp);
	}      
      d /= u.n();
      
      u.CompAdd(comp,-d);
    }
}

/*-----------------------------------------*/

