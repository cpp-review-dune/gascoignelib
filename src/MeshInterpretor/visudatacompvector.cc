#include  "visudatacompvector.h"

using namespace std;

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector()  
  : VisuData(), _v(NULL)
{
}

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector(const GlobalVector& v)  
  : VisuData(), _v(NULL)
{
  SetGlobalVector(&v);
}

/*---------------------------------------------------*/

void VisuDataCompVector::SetGlobalVector(const GlobalVector* v)
{
  assert(v);
  _v = v;
}

/*---------------------------------------------------*/

int    VisuDataCompVector::visucomp()     const 
{
  return _v->ncomp();
}

/*---------------------------------------------------*/

int    VisuDataCompVector::visun()        const 
{
  return _v->n();
}

/*---------------------------------------------------*/

double VisuDataCompVector::visudata(int i, int c) const 
{
  return (*_v)(i,c);
}
