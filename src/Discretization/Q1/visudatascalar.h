#ifndef __visudatascalar_h
#define __visudatascalar_h

#include  "visudata.h"
#include  "nvector.h"

/***************************************************************/

namespace Gascoigne
{
class VisuDataScalar : public VisuData
{
 protected:

  const DoubleVector& vR;

 public:

  virtual ~VisuDataScalar(){}
  VisuDataScalar(const DoubleVector& v) : vR(v) {}

  virtual int    visucomp()     const {return 1;}
  virtual int    visun()        const {return vR.size();}
  virtual double visudata(int i,int c) const { return vR[i];}
};
}

/***************************************************************/

#endif
