#include  "lpsintegrator.h"
#include  "lpsequation.h"
#include  "integrationformula.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

template<int DIM>
void LpsIntegrator<DIM>::Projection(const FemInterface& FEM) const
{
  if (DIM==2)
    {
      for(int ii=0; ii<FEM.n(); ii++) 
	{
	  FEM.init_test_functions(NLPS[ii],1.,ii);
	  // 	  NLPS[ii].m() = FEM.N  (ii);
// 	  NLPS[ii].x() = FEM.N_x(ii);
// 	  NLPS[ii].y() = FEM.N_y(ii);
	}
      NLPS[0].equ(-0.25, NLPS[4]);
      NLPS[2].equ(-0.25, NLPS[4]);
      NLPS[6].equ(-0.25, NLPS[4]);
      NLPS[8].equ(-0.25, NLPS[4]);
      
      NLPS[0].add(-0.5, NLPS[1], -0.5, NLPS[3]);
      NLPS[2].add(-0.5, NLPS[1], -0.5, NLPS[5]);
      NLPS[6].add(-0.5, NLPS[3], -0.5, NLPS[7]);
      NLPS[8].add(-0.5, NLPS[5], -0.5, NLPS[7]);
    }
  else if (DIM==3)
    {
      for(int ii=0; ii<FEM.n(); ii++) 
	{
	  FEM.init_test_functions(NLPS[ii],1.,ii);
// 	  NLPS[ii].m() = FEM.N  (ii);
// 	  NLPS[ii].x() = FEM.N_x(ii);
// 	  NLPS[ii].y() = FEM.N_y(ii);
// 	  NLPS[ii].z() = FEM.N_z(ii);
	}
      NLPS[0] .equ(-0.125, NLPS[13]);
      NLPS[2] .equ(-0.125, NLPS[13]);
      NLPS[6] .equ(-0.125, NLPS[13]);
      NLPS[8] .equ(-0.125, NLPS[13]);
      NLPS[18].equ(-0.125, NLPS[13]);
      NLPS[20].equ(-0.125, NLPS[13]);
      NLPS[24].equ(-0.125, NLPS[13]);
      NLPS[26].equ(-0.125, NLPS[13]);

      NLPS[0] .add(-0.5, NLPS[1] , -0.5, NLPS[3], -0.5, NLPS[9]);
      NLPS[2] .add(-0.5, NLPS[1] , -0.5, NLPS[5], -0.5, NLPS[11]);
      NLPS[6] .add(-0.5, NLPS[3] , -0.5, NLPS[7], -0.5, NLPS[15]);
      NLPS[8] .add(-0.5, NLPS[5] , -0.5, NLPS[7], -0.5, NLPS[17]);
      NLPS[18].add(-0.5, NLPS[19], -0.5, NLPS[21], -0.5, NLPS[9]);
      NLPS[20].add(-0.5, NLPS[19], -0.5, NLPS[23], -0.5, NLPS[11]);
      NLPS[24].add(-0.5, NLPS[21], -0.5, NLPS[25], -0.5, NLPS[15]);
      NLPS[26].add(-0.5, NLPS[23], -0.5, NLPS[25], -0.5, NLPS[17]);

      NLPS[0] .add(-0.25, NLPS[4] , -0.25, NLPS[10], -0.25, NLPS[12]);
      NLPS[2] .add(-0.25, NLPS[4] , -0.25, NLPS[10], -0.25, NLPS[14]);
      NLPS[6] .add(-0.25, NLPS[4] , -0.25, NLPS[16], -0.25, NLPS[12]);
      NLPS[8] .add(-0.25, NLPS[4] , -0.25, NLPS[16], -0.25, NLPS[14]);
      NLPS[18].add(-0.25, NLPS[22], -0.25, NLPS[10], -0.25, NLPS[12]);
      NLPS[20].add(-0.25, NLPS[22], -0.25, NLPS[10], -0.25, NLPS[14]);
      NLPS[24].add(-0.25, NLPS[22], -0.25, NLPS[16], -0.25, NLPS[12]);
      NLPS[26].add(-0.25, NLPS[22], -0.25, NLPS[16], -0.25, NLPS[14]);
    }
}

/*-----------------------------------------*/

template<int DIM>
void LpsIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  //assert(FEM.n()==9);
  NLPS.resize(FEM.n());
  MLPS.resize(FEM.n());

  const LpsEquation* LEQ = dynamic_cast<const LpsEquation*>(&EQ);
  assert(LEQ);

  const IntegrationFormulaInterface& IF = FormFormula();
  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      
      double vol = FEM.J();
      double weight  =  CellWeight * IF.w(k) * vol;
      FEM.x(x);
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      double h  = Volume2MeshSize(vol);
      LEQ->lpspoint(h,UH,QH,x);

      Projection(FEM);
      BasicIntegrator::universal_point(UHP,U,NLPS);
      for(int i=0;i<FEM.n();i++) MLPS[i].equ(weight,NLPS[i]);
            
      for (int i=0;i<FEM.n();i++)
	{
	  LEQ->StabForm(F.start(i),UH,UHP,MLPS[i]);
	}
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void LpsIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const
{
  assert(E.Ndof()==FEM.n());
  assert(E.Mdof()==FEM.n());
  assert(E.Ncomp()==U.ncomp());
  //assert(FEM.n()==9);
  NLPS.resize(FEM.n());
  MLPS.resize(FEM.n());

  const LpsEquation* LEQ = dynamic_cast<const LpsEquation*>(&EQ);
  assert(LEQ);

  const IntegrationFormulaInterface& IF = FormFormula();
  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      
      double vol = FEM.J();
      double weight  =  CellWeight * IF.w(k) * vol;
      FEM.x(x);
      BasicIntegrator::universal_point(FEM,UH,U);
      BasicIntegrator::universal_point(FEM,QH,Q);
      double h  = Volume2MeshSize(vol);
      LEQ->lpspointmatrix(h,UH,QH,x);

      Projection(FEM);

      for(int i=0;i<FEM.n();i++) MLPS[i].equ(weight,NLPS[i]);
      
      for (int j=0; j<FEM.n(); j++)
	{
	  for (int i=0; i<FEM.n(); i++)
	    {
	      E.SetDofIndex(i,j);
	      LEQ->StabMatrix(E, UH, NLPS[i], MLPS[j]);
	    }
	}
    }
}

/*-----------------------------------------------------------*/

template class LpsIntegrator<2>;
template class LpsIntegrator<3>;
}
