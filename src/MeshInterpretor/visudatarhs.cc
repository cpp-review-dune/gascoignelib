#include  "visudatarhs.h"


using namespace std;

/*-----------------------------------------*/

VisuDataRhs::VisuDataRhs(const RightHandSideData* rhs, const MeshInterface* m) 
  : RHS(rhs), M(m)
{
}

/*-----------------------------------------*/

int VisuDataRhs::visucomp() const
{
  return 1;
}

/*-----------------------------------------*/

int VisuDataRhs::visun() const
{
  return M->nnodes();
}

/*-----------------------------------------*/

double VisuDataRhs::visudata(int i,int c, const Vertex2d& v) const
{
  return (*RHS)(c,v);
}

/*-----------------------------------------*/

double VisuDataRhs::visudata(int i,int c, const Vertex3d& v) const
{
  return (*RHS)(c,v);
}
