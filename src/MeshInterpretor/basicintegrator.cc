#include  "basicintegrator.h"


using namespace std;

/*---------------------------------------------------------*/

BasicIntegrator::BasicIntegrator() : IntegratorInterface() 
{
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(const FemInterface& FEM, FemData& QH, const LocalData& Q) const
{
  QH.resize(Q.size());
  for(int i=0;i<Q.size();i++)
    {
      universal_point(FEM, QH[i], Q[i]);
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
