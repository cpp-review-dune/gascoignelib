#include  "localterminalcondition.h"


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

void LocalTerminalCondition::SetFemData(FemData& Q) const
{
  assert(Q.count("U"));
  q = &Q["U"];
}

/* ----------------------------------------- */

double LocalTerminalCondition::operator()(int c, const Vertex2d& v) const
{
  double x = v.x();
  double y = v.y();
  double eps = 1.e-3;

  if(c==0)
    {
      double uT = (*q)[0].m();

      x -= 100;
      y -= 150;
      double dist = (x*x+y*y)/(800*800+150*150);

      double ubar = pow((1.-dist),30.);
      ubar -= uT;

      return  ubar;
    }
  assert(0);
}
