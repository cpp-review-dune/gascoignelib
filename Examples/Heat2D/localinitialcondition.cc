#include  "localinitialcondition.h"


using namespace std;

/* ----------------------------------------- */

LocalInitialCondition::LocalInitialCondition(const LocalEquation* EQ)
{
  _us = EQ->GetUs();
  _vs = EQ->GetVs();
}

/* ----------------------------------------- */

double LocalInitialCondition::operator()(int c, const Vertex2d& v) const
{
  double x = v.x();
  double y = v.y();
  double eps1 = 2.e-7;
  double eps2 = 3.e-5;
  double eps3 = 1.2e-4;
  if(c==0)
    {
      double dist = - eps1*(x-0.1*y-225)*(x-0.1*y-675);
      return _us + dist;
    }
  else if(c==1)
    {
      return _vs - eps2*(x-450)-eps3*(y-150);
    }
  assert(0);
}
