#include "energyestimatorintegrator.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
void EnergyEstimatorIntegrator<DIM>::BasicInit()
{
  if (DIM==2)
  {
    IF = new QuadGauss4;
    assert(IF);
    _xi[0].x() = 0.;
    _xi[1].x() = 1.;
  }
  else
  {
    IF = new HexGauss8;
    assert(IF);
    _xi[0].x() = 0.; _xi[0].y() = 0.;
    _xi[1].x() = 1.; _xi[1].y() = 0.;
    _xi[2].x() = 1.; _xi[2].y() = 1.;
    _xi[3].x() = 0.; _xi[3].y() = 1.;
  }
}

/**********************************************************/

template<int DIM>
EnergyEstimatorIntegrator<DIM>::~EnergyEstimatorIntegrator()
{
  if(IF)
  {
    delete IF;
    IF = NULL;
  }
}

/**********************************************************/

template<int DIM>
void EnergyEstimatorIntegrator<DIM>::Jumps(LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile) const
{
  Vertex<DIM> n;

  F.ReInit(U.ncomp(),2*DIM-2);
  F.zero();

  for (int i=0; i<2*DIM-2; i++)
  {
    FEM.point_boundary(ile,_xi[i]);
    FEM.normal(n);

    BasicIntegrator::universal_point(FEM,UH,U);

    for (int c=0; c<U.ncomp(); c++)
    {
      F(i,c) += n.x() * UH[c].x() + n.y() * UH[c].y() + n.z() * UH[c].z();
    }
  }
}

/**********************************************************/

template<int DIM>
double EnergyEstimatorIntegrator<DIM>::JumpNorm(const FemInterface& FEM, fixarray<2*DIM-2,double> jumps, int ile) const
{
  double norm = 0.;
  for (int k=0; k<2*DIM-2; k++)
  {
    FEM.point_boundary(ile,_xi[k]);
    double h = Volume2MeshSize(FEM.J());
    double weight = (1.-DIM/4.)*h*FEM.G();
    norm += weight * jumps[k];
  }
  return norm;
}

/**********************************************************/

template<int DIM>
double EnergyEstimatorIntegrator<DIM>::Residual(const LocalVector& U, const FemInterface& FEM, const Equation& EQ, const DomainRightHandSide& RHS, const LocalNodeData& Q) const
{
  double res = 0.;
  DoubleVector F(U.ncomp());
  
  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<GetFormula().n(); k++)
  {
    GetFormula().xi(xi,k);
    FEM.point(xi);
    double vol = FEM.J();
    double h  = Volume2MeshSize(vol);
    double weight  = GetFormula().w(k) * vol;
    BasicIntegrator::universal_point(FEM,UH,U);
    BasicIntegrator::universal_point(FEM,QH,Q);
    FEM.x(x);
    RHS.SetFemData(QH);
    EQ.point(h,UH,QH,x);
    EQ.OperatorStrong(F,UH);
    double value = 0.;
    for (int c=0; c<U.ncomp(); c++)
    {
      value += (RHS(c,x)-F[c])*(RHS(c,x)-F[c]);
    }
    res += h*h*weight * value;
  }
  return res;
}

/**********************************************************/

template<int DIM>
double EnergyEstimatorIntegrator<DIM>::ResidualZeroRhs(const LocalVector& U, const FemInterface& FEM, const Equation& EQ, const LocalNodeData& Q) const
{
  double res = 0.;
  DoubleVector F(U.ncomp());

  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<GetFormula().n(); k++)
  {
    GetFormula().xi(xi,k);
    FEM.point(xi);
    double vol = FEM.J();
    double h  = Volume2MeshSize(vol);
    double weight  = GetFormula().w(k) * vol;
    BasicIntegrator::universal_point(FEM,UH,U);
    BasicIntegrator::universal_point(FEM,QH,Q);
    FEM.x(x);
    EQ.point(h,UH,QH,x);
    EQ.OperatorStrong(F,UH);
    double value = 0.;
    for (int c=0; c<U.ncomp(); c++)
    {
      value += F[c] * F[c];
    }
    res += h*h*weight * value;
  }
  return res;
}

/**********************************************************/

template class EnergyEstimatorIntegrator<2>;
template class EnergyEstimatorIntegrator<3>;
}
