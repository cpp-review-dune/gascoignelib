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
      xi[0].x() = 0.;
      xi[1].x() = 1.;
    }
  else if (DIM==3)
    {
      FormFormulaPointer() = new HexGauss8;
      ErrorFormulaPointer() = new HexGauss27;
      BoundaryFormulaPointer() = new QuadGauss4;
      xi[0].x() = 0.; xi[0].y() = 0.;
      xi[1].x() = 1.; xi[1].y() = 0.;
      xi[2].x() = 1.; xi[2].y() = 1.;
      xi[3].x() = 0.; xi[3].y() = 1.;
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
  F.ReInit(EQ.ncomp(),FEM.n());

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

template<int DIM>
double GalerkinIntegrator<DIM>::MeanMatrix(EntryMatrix& E, const FemInterface& FEM) const
{
  E.zero();

  const IntegrationFormulaInterface& IF = FormFormula();

  Vertex<DIM> xi;
  double omega = 0.;

  FemFunction NI(FEM.n()), NJ(FEM.n());
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      omega += weight;

      for(int i=0; i<FEM.n(); i++)
	{
	  NJ[i].m() = FEM.N(i);
	  NI[i].m() = weight * NJ[i].m();
	}
      for(int i=0; i<FEM.n(); i++)
	{
	  for(int j=0; j<FEM.n(); j++)
	    {
	      E.SetDofIndex(i,j);
	      E(0,0) += NI[i].m()*NJ[j].m();
	    }
	}
    }
  return omega;




//   int nv = MP->nodes_per_element(); // = 4 oder 8
//   std::vector<DerivativeVector> NI(nv), NJ(nv);
  
//   double omega = 0.;
//   Vertex2d xi;
//   for(int k=0;k<IF.n();++k)
//     {
//       IF.xi(xi,k);  // nur in 2d !!!!!!!!!!!!!!
//       FE.point(xi);
//       double  J = FE.J();
//       double  w =IF.w(k)*J;
//       omega += w;
      
//       for(int i=0;i<nv;i++)
// 	{
// 	  NJ[i].m() = FE.N  (i);
// 	  NI[i].m() = w * NJ[i].m();
// 	}
//       for(int i=0;i<nv;i++)
// 	{
// 	  for(int j=0;j<nv;j++)
// 	    {
// 	      E.SetDofIndex(i,j);
// 	      E(0,0) += NI[i].m()*NJ[j].m();
// 	    }
// 	}
//     }
//   return omega;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Jumps(LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile) const
{
  Vertex<DIM> n;

  F.ReInit(U.ncomp(),2*DIM-2);
  F.zero();

  for (int i=0; i<2*DIM-2; i++)
    {
      FEM.point_boundary(ile,xi[i]);
      FEM.normal(n);

      BasicIntegrator::universal_point(FEM,UH,U);

      for (int c=0; c<U.ncomp(); c++)
	{
	  F(i,c) += n.x() * UH[c].x() + n.y() * UH[c].y() + n.z() * UH[c].z();
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::JumpNorm(double& norm, const FemInterface& FEM, fixarray<2*DIM-2,double> jumps, int ile) const
{
  for (int k=0; k<2*DIM-2; k++)
    {
      FEM.point_boundary(ile,xi[k]);
      double h = Volume2MeshSize(FEM.J());
      double weight = (1.-DIM/4.)*h*FEM.G();
      norm += weight * jumps[k];
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Residual(double& res, const LocalVector& U, const FemInterface& FEM, const Equation& EQ, const RightHandSideData& RHS, const LocalData& Q) const
{
  DoubleVector F(U.ncomp());
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
      EQ.OperatorStrong(F,UH);
      double value = 0.;
      for (int c=0; c<U.ncomp(); c++)
	{
	  value += (RHS(c,x)-F[c])*(RHS(c,x)-F[c]);
	}
      res += h*h*weight * value;
    }
}

/* ----------------------------------------- */
/* ----------------------------------------- */

template class GalerkinIntegrator<2>;
template class GalerkinIntegrator<3>;
