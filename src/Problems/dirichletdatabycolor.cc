#include  "dirichletdatabycolor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
DirichletDataByColor::DirichletDataByColor(int c, const set<int>& cl, double s)
  : DirichletData(), comp(c), col(cl), scale(s)
{
}

/*-----------------------------------------*/

DirichletDataByColor::DirichletDataByColor(const vector<string>& args)
{
  int n = args.size();
  comp  = atoi(args[0].c_str());
  col.clear();
  for(int i=1;i<n-1;i++)  col.insert(atoi(args[i].c_str()));
  scale = atof(args[n-1].c_str());
}

/*-----------------------------------------*/

void DirichletDataByColor::operator()
(DoubleVector& b, const Vertex2d& v, int color) const
{
  b.zero();

  if(col.find(color)!=col.end()) b[comp] = scale;
}

/*-----------------------------------------*/

void DirichletDataByColor::operator()
(DoubleVector& b, const Vertex3d& v, int color) const
{
  b.zero();
  if(col.find(color)!=col.end()) b[comp] = scale;
}

/*-----------------------------------------*/

set<int> DirichletDataByColor::preferred_colors()const 
{
  return col;
}
}
