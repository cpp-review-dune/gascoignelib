#include  "galerkinintegrator.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
template<int DIM>
GalerkinIntegrator<DIM>::GalerkinIntegrator<DIM>() : BasicIntegrator(),
  IFF(0), IFE(0), IFB(0), IFM(0)
{
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new QuadGauss4;
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new QuadGauss9;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new LineGauss2;
      if (!MassFormulaPointer())     MassFormulaPointer() = new QuadGauss4;;
//      if (!MassFormulaPointer())     MassFormulaPointer() = new QuadTrapez;
    }
  else if (DIM==3)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new HexGauss8;
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new HexGauss27;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new QuadGauss4;
      if (!MassFormulaPointer())     MassFormulaPointer() = new HexGauss8;
//      if (!MassFormulaPointer())     MassFormulaPointer() = new HexTrapez;
    }
  assert(FormFormulaPointer());
  assert(ErrorFormulaPointer());
  assert(BoundaryFormulaPointer());
  assert(MassFormulaPointer());
}

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegrator<DIM>::~GalerkinIntegrator<DIM>()
{
  if(FormFormulaPointer())
  {
    delete FormFormulaPointer();
  }
  if(ErrorFormulaPointer())
  {
    delete ErrorFormulaPointer();
  }
  if(BoundaryFormulaPointer())
  {
    delete BoundaryFormulaPointer();
  }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Rhs(const DomainRightHandSide& f, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const
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
void GalerkinIntegrator<DIM>::BoundaryRhs(const BoundaryRightHandSide& f, LocalVector& F, const FemInterface& FEM, int ile, int col, const LocalNodeData& Q) const
{
  F.ReInit(f.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = BoundaryFormula();
  F.zero();
  Vertex<DIM> x, n;
  Vertex<DIM-1> xi;

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
void GalerkinIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  F.ReInit(EQ.GetNcomp(),FEM.n());

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
      EQ.SetFemData(QH);
      EQ.point(h,UH,x);
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NN,weight,i);
	  EQ.Form(F.start(i),UH,NN);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& Z, const LocalNodeData& Q) const
{
  F.ReInit(EQ.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  NNN.resize(FEM.n());
  EntryMatrix E;
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(Z.ncomp(),Z.ncomp());
  E.resize();
  E.zero();

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,UH,Z);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      // EQ.pointmatrix(h,QH["u"],QH,x);
      EQ.pointmatrix(h,QH["u"],x);
      double sw = sqrt(weight);
      for (int i=0; i<FEM.n(); i++)
	{
	  FEM.init_test_functions(NNN[i],sw,i);
	}
      for (int j=0; j<FEM.n(); j++)
	{
	  for (int i=0; i<FEM.n(); i++)
	    {
	      E.SetDofIndex(j,i);
	      EQ.Matrix(E,QH["u"],NNN[i],NNN[j]);
	    }
	}
    }
  for (int i=0; i<FEM.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  double sum = 0.;
	  for (int j=0; j<FEM.n(); j++)
	    {
	      for (int d=0; d<Z.ncomp(); d++)
		{
		  sum += E(j,i,d,c)*Z(j,d);
		}
	    }
	  F(i,c) += sum;
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile, int col, LocalNodeData& Q) const
{
  F.ReInit(BE.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = BoundaryFormula();

  F.zero();
  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      BE.SetFemData(QH);
      BE.pointboundary(h,UH,x,n);
      for (int i=0;i<FEM.n();i++)
        {
          FEM.init_test_functions(NN,weight,i);
          BE.Form(F.start(i),UH,NN,col);
        }
    }
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::MassMatrix(EntryMatrix& E, const FemInterface& FEM) const
{
  NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  int ncomp=1;
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();
  const IntegrationFormulaInterface& IF = MassFormula();

  Vertex<DIM> x, xi;
  double omega = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      omega += weight;
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NNN[i],1.,i);
	}
      for (int i=0;i<FEM.n();i++)
	{
	  for (int j=0;j<FEM.n();j++)
	    {
	      E.SetDofIndex(i,j);
	      E(0,0) += weight * NNN[j].m()*NNN[i].m();
	    }
	}
    }
  return omega;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
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
      EQ.SetFemData(QH);
      EQ.pointmatrix(h,UH,x);

      double sw = sqrt(weight);
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(NNN[i],sw,i);
	}
      for (int j=0; j<FEM.n(); j++)
	{
	  for (int i=0; i<FEM.n(); i++)
	    {
	      E.SetDofIndex(i,j);
	      EQ.Matrix(E,UH,NNN[j],NNN[i]);
	    }
	}
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void GalerkinIntegrator<DIM>::BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, int ile, int col, const LocalNodeData& Q) const
{
  NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(U.ncomp(),U.ncomp());
  E.resize();
  E.zero();

  const IntegrationFormulaInterface& IF = BoundaryFormula();

  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      BE.SetFemData(QH);
      BE.pointmatrixboundary(h,UH,x,n);
      double sw = sqrt(weight);
      for (int i=0;i<FEM.n();i++)
        {
          FEM.init_test_functions(NNN[i],sw,i);
        }
      for (int j=0;j<FEM.n();j++)
        {
          for (int i=0;i<FEM.n();i++)
            {
              E.SetDofIndex(i,j);
              BE.Matrix(E,UH,NNN[j],NNN[i],col);
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
void GalerkinIntegrator<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalNodeData& Q) const
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

/* ----------------------------------------- */
template<int DIM>
double GalerkinIntegrator<DIM>::ComputePointValue(const FemInterface& E, const Vertex<DIM>& p, const LocalVector& U, int comp) const
{
  E.point(p);
  BasicIntegrator::universal_point(E,UH,U);

  return UH[comp].m();
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  const IntegrationFormulaInterface& IF = FormFormula();

  Vertex<DIM> x, xi;
  double j = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      FEM.x(x);
      F.SetFemData(QH);
      j += weight * F.J(UH,x);
    }
  return j;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, const LocalNodeData& Q) const
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
	  // L2 Norm
	  dst(0,c) += weight * UH[c].m() * UH[c].m();
	  // H1 Seminorm
	  double a = UH[c].x()*UH[c].x() + UH[c].y()*UH[c].y();
	  if (DIM==3) a += UH[c].z() * UH[c].z();
	  dst(1,c) += weight * a;
	  // L8 Norm
	  dst(2,c) = Gascoigne::max(dst(2,c),fabs(UH[c].m()));
	}
    }
}

/* ----------------------------------------- */

template class GalerkinIntegrator<2>;
template class GalerkinIntegrator<3>;
}
