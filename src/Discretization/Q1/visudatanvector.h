#ifndef  __VisuDataNVector_h
#define  __VisuDataNVector_h


#include  "visudata.h"
#include  "nvector.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class VisuDataNVector : public VisuData
{
protected:

  const DoubleVector* vp;

public:


  VisuDataNVector(const DoubleVector& v) : vp(&v) {} 

  int    visucomp()     const {return 1;}
  int    visun()        const {return vp->size();}
  double visudata(int i,int c) const { return (*vp)[i];}

};
}

#endif
