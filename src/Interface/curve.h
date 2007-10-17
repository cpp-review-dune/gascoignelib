#ifndef  __Curve_h
#define  __Curve_h

#include "vertex.h"
#include "nvector.h"


namespace Gascoigne
{
/*-----------------------------------------*/


class Curve
{
 private:

 protected:

  
 public:

  Curve() {}
  virtual ~Curve() {}

  virtual int GetNcomp() const=0;

  virtual double X(double t) const = 0;
  virtual double Y(double t) const = 0;

  virtual double DX(double t) const = 0;
  virtual double DY(double t) const = 0;
  
  virtual Vertex2d operator()(double t) const
    {
      Vertex2d x;
      x.x()=X(t);
      x.y()=Y(t);
      return x;
    }
  
  virtual double NormD(double t) const {return sqrt(DX(t)*DX(t) + DY(t)*DY(t));}

  virtual void Normal(nvector<double>& n,double t) const
    {
      n.resize(2);
      double d = NormD(t);
      assert(d>=0);
      n[0] =  DY(t)/d;
      n[0] = -DX(t)/d;
    }

  virtual void Vertices(nvector<Vertex2d>& V,const nvector<double>& t) const
    {
      assert(V.size()==t.size());
      for(int i=0;i<t.size();++i)
	{
	  V[i].x() = X(t[i]);
	  V[i].y() = Y(t[i]);
	}
    }
};

}

#endif
