#include  "dirichletdatabycolor.h"

/*-----------------------------------------*/

DirichletDataByColor::DirichletDataByColor(int c, const std::set<int>& cl, double s)
  : DirichletData(), comp(c), col(cl), scale(s)
{
}

/*-----------------------------------------*/

DirichletDataByColor::DirichletDataByColor(const std::vector<std::string>& args)
{
  int n = args.size();
  comp  = atoi(args[0].c_str());
  col.clear();
  for(int i=1;i<n-1;i++)  col.insert(atoi(args[i].c_str()));
  scale = atof(args[n-1].c_str());
}

/*-----------------------------------------*/

void DirichletDataByColor::operator()
(Vector& b, const Vertex2d& v, int color) const
{
  b.zero();

  if(col.find(color)!=col.end()) b[comp] = scale;
}

/*-----------------------------------------*/

void DirichletDataByColor::operator()
(Vector& b, const Vertex3d& v, int color) const
{
  b.zero();
  if(col.find(color)!=col.end()) b[comp] = scale;
}

/*-----------------------------------------*/

std::set<int> DirichletDataByColor::preferred_colors()const 
{
  return col;
}
