#include  "benchmark3ddirichletdata.h"
#include  "usefullfunctionsbd.h"

using namespace std;

/* ----------------------------------------- */

BenchMark3dDirichletData::BenchMark3dDirichletData(const std::string& paramfile)
{
  vmax = 0.45;
}

/* ----------------------------------------- */

void BenchMark3dDirichletData::operator()(Vector& b, const Vertex3d& v,int color) const
{
  double x = v.x();  double y = v.y(); double z = v.z();

  b.zero();
  if (x==0.)
    {
      double high = 4.1;
      b[1] = vmax * ParabelFunction(y,0.,high) * ParabelFunction(z,0.,high);
    }
}
