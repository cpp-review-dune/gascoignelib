#include  "visudatacompvector.h"

using namespace std;

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector()  
: VisuData()
{
  isizes.push_back(0);
}

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector(const GlobalVector& v)  
: VisuData()
{
  isizes.push_back(0);
  AddGlobalVector(&v);
}

/*---------------------------------------------------*/

void VisuDataCompVector::AddGlobalVector(const GlobalVector* v)
{
  assert(v);

  vvp.push_back(v);
  assert(vvp.size());
  assert(vvp[0]);
  if(v->n()!=vvp[0]->n())
    {
      cerr << "VisuDataCompVector::AddVector()" << endl;
      cerr << "n unterschiedlich gross\n";
      abort();
    }
  int n = isizes[isizes.size()-1];
  isizes.push_back(n+v->ncomp());
//   cerr << "isizes\n";
//   copy(isizes.begin(),isizes.end(),ostream_iterator<int>(cerr," "));
//   cerr << endl;
}

/*---------------------------------------------------*/

int    VisuDataCompVector::visucomp()     const 
{
  return isizes[isizes.size()-1];
}

/*---------------------------------------------------*/

int    VisuDataCompVector::visun()        const 
{
  return vvp[0]->n();
}

/*---------------------------------------------------*/

pair<int,int> VisuDataCompVector::GetIndex(int c) const
{
  for(int i=1;i<isizes.size();i++)
    {
      if(c<isizes[i]) return make_pair(i-1,c-isizes[i-1]);
    }
  cerr << "Not Found\n";
  abort();
}

/*---------------------------------------------------*/

double VisuDataCompVector::visudata(int i,int c) const 
{ 
  pair<int,int> p = GetIndex(c);
  return *(vvp[p.first]->start(i) + p.second);
}
