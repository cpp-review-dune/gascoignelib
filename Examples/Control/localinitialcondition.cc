#include  "localinitialcondition.h"


using namespace std;

/* ----------------------------------------- */

LocalInitialCondition::LocalInitialCondition(const LocalEquation* EQ)
{
  cerr << "equation is: " << EQ->GetName() << endl;
  ncomp = EQ->ncomp();
}

/* ----------------------------------------- */

double LocalInitialCondition::operator()(int c, const Vertex2d& v) const
{
  double x = v.x();
  double y = v.y();

  if(c==0)
    {
      x -= 800;
      y -= 150;
      double dist = (x*x+y*y)/(800*800+150*150);
      return pow((1.-dist),30.);
    }
  assert(0);
}
