#ifndef  __ConstantRightHandSide_h
#define  __ConstantRightHandSide_h

#include  "righthandsidedata.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class OneRightHandSideData : public RightHandSideData
{
  protected :

    int ncomp;

  public  : 

    OneRightHandSideData(int n) : ncomp(n), RightHandSideData() {}
    std::string GetName() const { return "one" ;} 
    int GetNcomp() const {return ncomp;}
    double operator()(int c, const Vertex2d& v)const {return 1.;}
    double operator()(int c, const Vertex3d& v)const {return 1.;}
};

/*-----------------------------------------*/

class ConstantRightHandSideData : public RightHandSideData
{
protected :
  int     comp,ncomp;
  double  d;

public  : 
  ConstantRightHandSideData(const std::vector<std::string>& args);
  std::string GetName() const {return "constant";} 
  int GetNcomp() const {return ncomp;}
  double operator()(int c, const Vertex2d& v)const;
  double operator()(int c, const Vertex3d& v)const;
};

/*-----------------------------------------*/

class OneComponentRightHandSideData : public RightHandSideData
{
 protected:
  
  int  ncomp, comp;  // ist die Komponente die Eins ist

 public:
  
  OneComponentRightHandSideData(int n, int c) : 
    RightHandSideData(), ncomp(n), comp(c) {}

  std::string GetName() const {return "one_onecomp";} 

  int GetNcomp() const { return ncomp;}

  double operator()(int c, const Vertex2d&)const 
    {
      if (c==comp) return 1.;
      return 0.;
    }
  double operator()(int c, const Vertex3d&)const 
    {
      if (c==comp) return 1.;
      return 0.;
    }
};

/*-----------------------------------------*/

class RectangleRightHandSideData : public RightHandSideData
{
  int      ncomp, comp;  // ist die Komponente die Eins ist
  double   x0, x1, y0, y1, z0, z1;

public:

  RectangleRightHandSideData(int n, int c, 
			     double xx0, double xx1, double yy0, double yy1) : 
    RightHandSideData(), ncomp(n), comp(c) 
    { 
      x0 = xx0; x1 = xx1; y0 = yy0; y1 = yy1;
      z0 = z1 = 0.;
    }
  RectangleRightHandSideData(int n, int c, 
			     double xx0, double xx1, double yy0, double yy1,
			     double zz0, double zz1) : 
    RightHandSideData(), ncomp(n), comp(c) 
    { 
      x0 = xx0; x1 = xx1; y0 = yy0; y1 = yy1;
      z0 = zz0; z1 = zz1;
    }

  std::string GetName() const {return "RectangleRightHandSideData";} 

  int GetNcomp() const {return ncomp;}

  double operator()(int c, const Vertex2d& V)const 
    {
      if (c!=comp) return 0.;
      if ((V.x()>x1) || (V.x()<x0)) return 0.;
      if ((V.y()>y1) || (V.y()<y0)) return 0.;
      return 1.;
    }
  double operator()(int c, const Vertex3d& V)const 
    {
      if (c!=comp) return 0.;
      if ((V.x()>x1) || (V.x()<x0)) return 0.;
      if ((V.y()>y1) || (V.y()<y0)) return 0.;
      if ((V.z()>z1) || (V.z()<z0)) return 0.;
      return 1.;
    }
};
}

#endif
