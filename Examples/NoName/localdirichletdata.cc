#include  "localdirichletdata.h"
#include  "usefullfunctionsbd.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

LocalDirichletData::LocalDirichletData(const ParamFile* paramfile)
{
  vmax = 0.3;
}

/* ----------------------------------------- */

void LocalDirichletData::operator()(Vector& b, const Vertex2d& v,int color) const
{
  b.zero();
  if (color!=4)
    {
      b[1] = vmax;
    }
}

/* ----------------------------------------- */
