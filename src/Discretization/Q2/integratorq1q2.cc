#include "integratorq1q2.h"
#include "patchintegrationformula.h"

namespace Gascoigne
{
/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, const LocalNodeData& Q) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(U.ncomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;

  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point(xi);
      FemL.point(xi);
      double vol = FemL.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemL,UH,U);
      BasicIntegrator::universal_point(FemL,QH,Q);
      FemL.x(x);
      EQ.SetFemData(QH);
      EQ.point(h,UH,x);
      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(NN,weight,i);
	  EQ.Form(F.start(i),UH,NN);
	}
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalNodeData& Q) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(RHS.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;
  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point(xi);
      FemL.point(xi);
      double vol = FemL.J();
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemL,QH,Q);
      RHS.SetFemData(QH);
      FemL.x(x);

      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(NN,weight,i);
	  RHS(F.start(i),NN,x);
	}
    } 
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalNodeData& Q) const
{
  b.zero();

  Vertex<DIM> x;
  E.point(p);     
  E.x(x);
  BasicIntegrator::universal_point(E,QH,Q);
  DRHS.SetFemData(QH);

  for (int i=0; i<E.n(); i++)
    {
      E.init_test_functions(NN,1.,i);
      DRHS.operator()(j,b.start(i),NN,x);
    }
}

/*---------------------------------------------------*/

template IntegratorQ1Q2<2>;
template IntegratorQ1Q2<3>;
}
