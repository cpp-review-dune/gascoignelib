#ifndef __integrationformulainterface_h
#define __integrationformulainterface_h

#include  "nvector.h"
#include  "vertex.h"

/*------------------------------------------------------------*/

class IntegrationFormulaInterface
{
public:

  IntegrationFormulaInterface() {}

  virtual ~IntegrationFormulaInterface() {}

  virtual int    n()      const=0;
  virtual double w(int k) const=0;

  virtual void xi(Vertex1d& v, int k) const {assert(0);}
  virtual void xi(Vertex2d& v, int k) const {assert(0);}
  virtual void xi(Vertex3d& v, int k) const {assert(0);}
};

/*------------------------------------------------------------*/

#endif
