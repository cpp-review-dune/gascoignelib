#include  "benchmarkdirichletdata.h"
#include  "usefullfunctionsbd.h"

using namespace std;

/* ----------------------------------------- */

BenchMarkDirichletData::BenchMarkDirichletData(const std::string& paramfile)
{
  vmax = 0.3;
}

/* ----------------------------------------- */

void BenchMarkDirichletData::operator()(Vector& b, const Vertex2d& v,int color) const
{
  double x = v.x();  double y = v.y();

  b.zero();
  if (color!=80)
    {
      double high = 4.1;
      b[1] = vmax * ParabelFunction(y,0.,high);
    }
  if(b.size()>3)
    {
      double tmin = 300.;
      double tmax = 600.;
      //double tmax = 300.;
      
      double x0 = 2.;
      double y0 = 2.;
      
      if ((x-x0)*(x-x0)+(y-y0)*(y-y0) > 1. ) b[3]= tmin;
      else b[3] = tmax;
    }
  if(b.size()>4)
    {
      double T = b[3];
      b[4] = 273./ T;
    }
}

/* ----------------------------------------- */
