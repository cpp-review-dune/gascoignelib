#include  "visudatacompvector.h"

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector()  
: VisuData()
{
  isizes.push_back(0);
}

/*---------------------------------------------------*/

VisuDataCompVector::VisuDataCompVector(const CompVector<double>& v)  
: VisuData()
{
  isizes.push_back(0);
  AddVector(v);
}

/*---------------------------------------------------*/

void VisuDataCompVector::AddVector(const CompVector<double>& v)
{
  vvp.push_back(&v);
  assert(vvp.size());
  assert(vvp[0]);
  if(v.n()!=vvp[0]->n())
    {
      std::cerr << "VisuDataCompVector::AddVector()" << std::endl;
      std::cerr << "n unterschiedlich gross\n";
      abort();
    }
  int n = isizes[isizes.size()-1];
  isizes.push_back(n+v.ncomp());
//   std::cerr << "isizes\n";
//   copy(isizes.begin(),isizes.end(),std::ostream_iterator<int>(std::cerr," "));
//   std::cerr << std::endl;
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

std::pair<int,int> VisuDataCompVector::GetIndex(int c) const
{
  for(int i=1;i<isizes.size();i++)
    {
      if(c<isizes[i]) return std::make_pair(i-1,c-isizes[i-1]);
    }
  std::cerr << "Not Found\n";
  abort();
}

/*---------------------------------------------------*/

double VisuDataCompVector::visudata(int i,int c) const 
{ 
  std::pair<int,int> p = GetIndex(c);
  return *(vvp[p.first]->start(i) + p.second);
}
