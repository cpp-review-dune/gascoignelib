#include  "galerkinintegrator.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegrator<DIM>::GalerkinIntegrator<DIM>() : BasicIntegrator() 
{
  if (DIM==2)
    {
      FormFormulaPointer() = new QuadGauss4;
      ErrorFormulaPointer() = new QuadGauss9;
      BoundaryFormulaPointer() = new LineGauss2;
    }
  else if (DIM==3)
    {
      FormFormulaPointer() = new HexGauss8;
      ErrorFormulaPointer() = new HexGauss27;
      BoundaryFormulaPointer() = new QuadGauss4;
    }
  assert(FormFormulaPointer());
  assert(ErrorFormulaPointer());
  assert(BoundaryFormulaPointer());
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Rhs(const RightHandSideData& f, LocalVector& F, const FemInterface& FEM, const LocalData& Q) const
{
  F.ReInit(f.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = FormFormula();

  F.zero();
  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,QH,Q);
      f.SetFemData(QH);
      FEM.x(x);
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NN,weight,i);
	  f(F.start(i),NN,x);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::RhsNeumann(const NeumannData& f, LocalVector& F, const FemInterface& FEM, int ile, int col, const LocalData& Q) const
{
  assert(DIM==2);
  F.ReInit(f.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = BoundaryFormula();
  F.zero();
  Vertex<DIM> x, n;
  Vertex1d xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,QH,Q);
      f.SetFemData(QH);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NN,weight,i);
	  f(F.start(i),NN,x,n,col);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const
{
  F.ReInit(U.ncomp(),FEM.n());

  const IntegrationFormulaInterface& IF = FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      EQ.point(h,UH,QH,x);
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NN,weight,i);
	  EQ.Form(F.start(i),UH,NN);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::MassMatrix(EntryMatrix& E, const FemInterface& FEM) const
{
  E.SetDimensionDof(FEM.n(),FEM.n());
  int ncomp=1;
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();

  const IntegrationFormulaInterface& IF = FormFormula();

  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      for (int i=0;i<FEM.n();i++)
	{
	  for (int j=0;j<FEM.n();j++)
	    {
	      E.SetDofIndex(i,j);
	      E(0,0) += weight*FEM.N(i)*FEM.N(j);
	    }
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const
{
  NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(U.ncomp(),U.ncomp());
  E.resize();
  E.zero();

  const IntegrationFormulaInterface& IF = FormFormula();

  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      EQ.pointmatrix(h,UH,QH,x);

      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NNN[i],weight,i);
	}
      for (int j=0;j<FEM.n();j++)
	{
	  FEM.init_test_functions(MM,1.,j);
	  for (int i=0;i<FEM.n();i++)
	    {
	      E.SetDofIndex(i,j);
	      EQ.Matrix(E,UH,MM,NNN[i]);
	    }
	}
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void GalerkinIntegrator<DIM>::RhsPoint
(LocalVector& F, const FemInterface& E, const Vertex<DIM>& p, int comp) const
{
  F.zero();

  E.point(p);
  for (int i=0; i<E.n(); i++)
    {
      E.init_test_functions(NN,1.,i);
      F(i,comp) += NN.m();
    }
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const
{
  const IntegrationFormulaInterface& IF = FormFormula();

  Vertex<DIM> x, xi;
  double j = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      j += weight * F.J(UH,x);
    }
  return j;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, const LocalData& Q) const
{
  const IntegrationFormulaInterface& IF = ErrorFormula();

  Vertex<DIM> x, xi;
  dst.zero();
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FE.point(xi);
      BasicIntegrator::universal_point(FE,UH,U);
      double vol = FE.J();
      double weight = IF.w(k) * vol;

      FE.x(x);
      for(int c=0; c<U.ncomp(); c++)
	{
	  UH[c].m() -= ES(c,x);
	  UH[c].x() -= ES.x(c,x);
	  UH[c].y() -= ES.y(c,x);
	  if (DIM==3) UH[c].z() -= ES.z(c,x);
	}
      for(int c=0; c<U.ncomp(); c++) 
	{
	  UH[c].m() *= weight*UH[c].m();
	  UH[c].x() *= weight*UH[c].x();
	  UH[c].y() *= weight*UH[c].y();
	  if (DIM==3) UH[c].z() *= weight*UH[c].z();
	}
      for(int c=0; c<U.ncomp(); c++) 
	{
	  dst(0,c) += UH[c].m();
	  dst(1,c) += UH[c].x() + UH[c].y();
	  if (DIM==3) dst(1,c) += UH[c].z();
	  dst(2,c) = GascoigneMath::max(dst(2,c),sqrt(UH[c].m()/weight));
	}
    }
}

/* ----------------------------------------- */

template GalerkinIntegrator<2>;
template GalerkinIntegrator<3>;
