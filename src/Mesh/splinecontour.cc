#include "splinecontour.h"

/*****************************************************************************/

SplineCoord::SplineCoord(int n) : 
  xa(n), ha(n) 
{}

/*****************************************************************************/

SplineContour::SplineContour(const SplineCoord* p) : spline(p->x(),p->y()),
  cp(p)
{
}

/*****************************************************************************/

SplineContour::~SplineContour()
{
  if (cp!=NULL) delete cp;
  cp = NULL;
}

/*****************************************************************************/

double SplineContour::operator()(const numfixarray<2,double>& c) const 
{
  double r = cp->right();
  double l = cp->left();
  if ( (c[0]>l) && (c[0]<r) )
    {
      return c[1] - spline(fabs(c[0]));
    }
  return 0.;
}

/*****************************************************************************/
