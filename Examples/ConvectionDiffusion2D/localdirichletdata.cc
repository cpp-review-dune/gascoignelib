#include  "localdirichletdata.h"


using namespace std;

/* ----------------------------------------- */

void LocalDirichletData::operator()(Vector& b, const Vertex2d& v, int col) const
{
  if(col!=80)
    {
      b[0] = 0.;
    }
  else
    {
      b[0] = 1.;
    }
}
