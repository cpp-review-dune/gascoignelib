#ifndef  __VisuDataNVector_h
#define  __VisuDataNVector_h


#include  "visudata.h"
#include  "nvector.h"

/*-----------------------------------------*/


class VisuDataNVector : public VisuData
{
protected:

  typedef  nvector<double>  Vector;

  const Vector* vp;

public:


  VisuDataNVector(const Vector& v) : vp(&v) {} 

  int    visucomp()     const {return 1;}
  int    visun()        const {return vp->size();}
  double visudata(int i,int c) const { return (*vp)[i];}

};


#endif
