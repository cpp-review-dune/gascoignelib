#include  "visudatacompvectorfunction.h"

using namespace std;

static MySolutionFunction MYSF;

/*---------------------------------------------------*/

double MySolutionFunction::operator()(int i, const CompVector<double>& U) const 
{
  double norm = U(i,1)*U(i,1) + U(i,2)*U(i,2);
  double T = U(i,3);
  //double R = 2.87e-6;
  double R = 1.646e-2;//7.316e-3;
  double gamma = 1.;
  return sqrt( norm / (gamma*R*T));
};

/*---------------------------------------------------*/

VisuDataCompVectorFunction::VisuDataCompVectorFunction()  
: VisuDataCompVector()
{
  outcome.resize(0);
  outcome.push_back(&MYSF);
}

/*---------------------------------------------------*/

VisuDataCompVectorFunction::VisuDataCompVectorFunction(const CompVector<double>& v)  
: VisuDataCompVector(v)
{
  outcome.resize(0);
  outcome.push_back(&MYSF);
}

/*---------------------------------------------------*/

int VisuDataCompVectorFunction::visucomp() const 
{
  return VisuDataCompVector::visucomp() + outcome.size();
}

/*---------------------------------------------------*/

pair<int,int> VisuDataCompVectorFunction::GetIndex(int c) const
{
  for(int i=0;i<isizes.size()-1;i++)
    {
      if(c<isizes[i+1]+1) return make_pair(i,c-isizes[i]);
    }
  make_pair(0,c+100);
}

/*---------------------------------------------------*/

double VisuDataCompVectorFunction::visudata(int i,int c) const 
{ 
  pair<int,int> p = GetIndex(c);

  int comp = p.second;


  if (comp<isizes[isizes.size()-1])
    {
      return *(vvp[p.first]->start(i) + comp);
    }
  else
    {
      return (*outcome[0])(i,*vvp[p.first]);
    }
}
