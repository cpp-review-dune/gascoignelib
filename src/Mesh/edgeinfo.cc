#include "edgeinfo.h"

using namespace std;
using namespace Gascoigne;

/**********************************************************/

template<int DIM>
void EdgeInfo<DIM>::basicInit(const Edge* edge, int ncomp, const fixarray<2*DIM-2,int>& vertex)
{
  _count  = 0;
  _edge   = edge;
  _vertex = vertex;
  _u.ReInit(ncomp,2*DIM-2);
  _u.zero();
}

/**********************************************************/

template<int DIM>
void EdgeInfo<DIM>::addNodes(const LocalVector& u)
{
  for (int i=0; i<2*DIM-2; i++)
    {
      for (int c=0; c<_u.ncomp(); c++)
	{
	  _u(i,c) += u(i,c);
	}
    }
  _count++;
}

/**********************************************************/

template<int DIM>
void EdgeInfo<DIM>::showStatistics() const
{
  cout << _count << ": " << _vertex << " = " << _u << endl;
}

/**********************************************************/

template<int DIM>
const LocalVector& EdgeInfo<DIM>::getValue() const
{
  return _u;
}

/**********************************************************/

template<int DIM>
fixarray<2*DIM-2,double> EdgeInfo<DIM>::getNorm() const
{
  fixarray<2*DIM-2,double> norm(0.);

  for (int i=0; i<2*DIM-2; i++)
    {
      for (int c=0; c<_u.ncomp(); c++)
	{
	  norm[i] += _u(i,c) * _u(i,c);
	}
    }
  return norm;
}

/**********************************************************/

template<int DIM>
const Edge& EdgeInfo<DIM>::getEdge() const
{
  return *_edge;
}

/**********************************************************/

template<int DIM>
const fixarray<2*DIM-2,int>& EdgeInfo<DIM>::vertex() const
{
  return _vertex;
}

/**********************************************************/

template<int DIM>
int EdgeInfo<DIM>::getCount() const
{
  return _count;
}

/**********************************************************/

template EdgeInfo<2>;
template EdgeInfo<3>;
