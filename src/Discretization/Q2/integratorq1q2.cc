#include "integratorq1q2.h"
#include "patchintegrationformula.h"

namespace Gascoigne
{
/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::BasicInit()
{
}

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
void IntegratorQ1Q2<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& Z, const LocalNodeData& Q) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(Z.ncomp(),FemH.n());
  F.zero();

  NNN.resize(FemH.n());
  EntryMatrix E;
  E.SetDimensionDof(FemH.n(),FemL.n());
  E.SetDimensionComp(Z.ncomp(),Z.ncomp());
  E.resize();
  E.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;
  TestFunction MM;
  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point(xi);
      FemL.point(xi);
      double vol = FemL.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemL,UH,Z);
      BasicIntegrator::universal_point(FemL,QH,Q);
      FemL.x(x);
      //EQ.pointmatrix(h,QH["u"],QH,x);
      EQ.pointmatrix(h,QH["u"],x);
      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(NNN[i],weight,i);
	}
      for (int j=0;j<FemL.n();j++)
	{
	  FemL.init_test_functions(MM,1.,j);
	  for (int i=0;i<FemH.n();i++)
	    {
	      E.SetDofIndex(j,i);
	      EQ.Matrix(E,QH["u"],NNN[i],MM);
	    }
	}
    }
  for (int i=0;i<FemH.n();i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  double sum = 0.;
	  for (int j=0; j<FemL.n(); j++)
	    {
	      for (int d=0; d<Z.ncomp(); d++)
		{
		  sum += E(j,i,d,c)*Z(j,d);
		}
	    }
	  F(i,c) += sum;
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
void IntegratorQ1Q2<DIM>::BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, 
    const FemInterface& FemL, int ile, int col, const LocalNodeData& Q) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(RHS.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula1d<3,LineGauss3>;
  else        IF = new PatchFormula2d<9,QuadGauss9>;

  Vertex<DIM> x, n;
  Vertex<DIM-1> xi;

  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point_boundary(ile,xi);
      FemL.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FemL,QH,Q);
      RHS.SetFemData(QH);
      FemL.x(x);
      FemL.normal(n);
      double  h = FemL.G();
      double  weight = IF->w(k)*h;
      for (int i=0;i<FemH.n();i++)
      {
        FemH.init_test_functions(NN,weight,i);
        RHS(F.start(i),NN,x,n,col);
      }
    }
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, 
    const FemInterface& FemL, const LocalVector& U, int ile, int col, LocalNodeData& Q) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(BE.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula1d<3,LineGauss3>;
  else        IF = new PatchFormula2d<9,QuadGauss9>;
  
  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;

  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point_boundary(ile,xi);
      FemL.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FemL,UH,U);
      BasicIntegrator::universal_point(FemL,QH,Q);
      FemL.x(x);
      FemL.normal(n);
      double  h = FemL.G();
      double  weight = IF->w(k)*h;
      BE.SetFemData(QH);
      BE.pointboundary(h,UH,x,n);
      for (int i=0;i<FemH.n();i++)
      {
        FemH.init_test_functions(NN,weight,i);
        BE.Form(F.start(i),UH,NN,col);
      }
    }
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
