#ifndef __integrationformulainterface_h
#define __integrationformulainterface_h

#include  "nvector.h"
#include  "vertex.h"

/*------------------------------------------------------------*/

class IntegrationFormulaInterface
{
protected:

  int              in;
  nvector<double>  iw;

public:

  IntegrationFormulaInterface() {}
  IntegrationFormulaInterface(int n) {}
  IntegrationFormulaInterface(const IntegrationFormulaInterface& IF) : iw(IF.w())
    {
      in = IF.n();
    }

  virtual ~IntegrationFormulaInterface() {}

  void init(int n);

  virtual int    n()                 const { return in;}
  virtual double w(int k)            const { return iw[k];}
  virtual const nvector<double>& w() const { return iw;}

  virtual void xi(Vertex1d& v, int k) const {}
  virtual void xi(Vertex2d& v, int k) const {}
  virtual void xi(Vertex3d& v, int k) const {}
};

/*------------------------------------------------------------*/

#endif
