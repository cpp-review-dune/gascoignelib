#ifndef __splinecontour_h
#define __splinecontour_h

#include  "boundaryfunction.h"
#include  "spline.h"

/*******************************************************/

class SplineCoord
{
 protected:

  nvector<double>   xa, ha;
  double            leftbound, rightbound;

 public:

  SplineCoord(int n);

  const nvector<double>& x() const { return xa; }
  const nvector<double>& y() const { return ha; }
  
  double left () const { return leftbound;}
  double right() const { return rightbound;}
};

/*******************************************************/

class SplineContour : public BoundaryFunction<2>
{
  const SplineCoord*  cp;
  CubicSpline         spline;
  
  public :
    
  SplineContour(const SplineCoord*);
  ~SplineContour();

  double operator()(const numfixarray<2,double>& c) const;
};

/********************************************************************/

#endif
