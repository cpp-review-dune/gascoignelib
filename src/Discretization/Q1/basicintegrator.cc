#include  "basicintegrator.h"


using namespace std;

/*---------------------------------------------------------*/

namespace Gascoigne
{
BasicIntegrator::BasicIntegrator() : IntegratorInterface() 
{
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(const FemInterface& FEM, FemData& QH, const LocalNodeData& Q) const
{
  QH.clear();
  LocalNodeData::const_iterator p=Q.begin();
  for(; p!=Q.end(); p++)
    {
      universal_point(FEM, QH[p->first], p->second);
    }
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(const FemInterface& FEM, FemFunction& UH, const LocalVector& U) const
{
  UH.resize(U.ncomp());

  for (int c=0; c<UH.size(); c++)  UH[c].zero();

  for (int i=0; i<FEM.n(); i++)
    {
      FEM.init_test_functions(NN,1.,i);
      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].add(U(i,c),NN);
	}
    }
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(FemFunction& UH, const LocalVector& U, const FemFunction& NN) const
{
  UH.resize(U.ncomp());
  for (int c=0; c<UH.size(); c++)  UH[c].zero();

  for (int i=0; i<NN.size(); i++)
    {
      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].add(U(i,c),NN[i]);
	}
    }
}
}
