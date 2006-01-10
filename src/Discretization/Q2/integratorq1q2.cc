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

/**
 * Wandelt die lokale Nummerierung der Zellen eines Patches von der Gascoigne-nummerierung in die 
 * vom Integrator verwendete um! Funktioniert aus Symmetrie-Gruenden auch andersherum.
 * Erwartet eine Zahl zwischen 0 und 3 (2-D) bzw 0 bis 7 (3-D)
*/

template<int DIM>
int IntegratorQ1Q2<DIM>::PatchMeshNr2IntegratorNr(int in) const
{
    if(DIM==2)
    {
	return (1-(in>>1))*in+(in>>1)*(3-(in%2));
    }
    if(DIM==3)
    {
	int tmp= in-((in>>2)<<2);
	return ((in>>2)<<2)+(1-(tmp>>1))*tmp+(tmp>>1)*(3-(tmp%2));
    }
    std::cerr<<"Integrator Q1/Q2: Dimension "<<DIM<<" nicht implementiert!"<<std::endl;
    abort();
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
    const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(U.ncomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;

  int numcells;
  if(DIM==2) numcells = 4;
  else numcells = 8;
  
  //Soviele Integr. Pkte je Zelle
  numcells = IF->n()/numcells;
  //Test, ist das sinnvollgewesen
  if(IF->n()%numcells != 0)
  {
      std::cerr<<"Integrator Q1/Q2: Falsche Anzahl Zellen je Patch oder unzulaessige Integrationsformel!"<<std::endl;
      abort();
  }

  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point(xi);
      FemL.point(xi);
      double vol = FemL.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemL,_UH,U);
      BasicIntegrator::universal_point(FemL,_QH,Q);
      FemL.x(x);

      EQ.SetFemData(_QH);
      //Die Zelldatensetzen
      if(k%numcells == 0)
      {
	  //Wir sind am anfang einer neuen Zelle
	  //Die CellData aus den LocalCellData hohlen.
	  int IntegrCellNum = k/numcells;
	  int PatchCellNum=PatchMeshNr2IntegratorNr(IntegrCellNum);
	  
	  universal_point(_QCH,QC,PatchCellNum);
	  EQ.SetCellData(_QCH);
      }
      
      EQ.point(h,_UH,x);

      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(_NN,weight,i);
	  EQ.Form(F.start(i),_UH,_NN);
	}
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
    const LocalVector& Z, const LocalNodeData& Q, const LocalCellData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(Z.ncomp(),FemH.n());
  F.zero();

  _NNN.resize(FemH.n());
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
      BasicIntegrator::universal_point(FemL,_UH,Z);
      BasicIntegrator::universal_point(FemL,_QH,Q);
      FemL.x(x);
      //EQ.pointmatrix(h,_QH["u"],_QH,x);
      EQ.pointmatrix(h,_QH["u"],x);
      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(_NNN[i],weight,i);
	}
      for (int j=0;j<FemL.n();j++)
	{
	  FemL.init_test_functions(MM,1.,j);
	  for (int i=0;i<FemH.n();i++)
	    {
	      E.SetDofIndex(j,i);
	      EQ.Matrix(E,_QH["u"],_NNN[i],MM);
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
void IntegratorQ1Q2<DIM>::Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, 
    const FemInterface& FemL, const LocalNodeData& Q, const LocalCellData& QC) const
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
      BasicIntegrator::universal_point(FemL,_QH,Q);
      RHS.SetFemData(_QH);
      FemL.x(x);

      for (int i=0;i<FemH.n();i++)
	{
	  FemH.init_test_functions(_NN,weight,i);
	  RHS(F.start(i),_NN,x);
	}
    } 
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, 
    const FemInterface& FemL, int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const
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
      BasicIntegrator::universal_point(FemL,_QH,Q);
      RHS.SetFemData(_QH);
      FemL.x(x);
      FemL.normal(n);
      double  h = FemL.G();
      double  weight = IF->w(k)*h;
      for (int i=0;i<FemH.n();i++)
      {
        FemH.init_test_functions(_NN,weight,i);
        RHS(F.start(i),_NN,x,n,col);
      }
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, 
    const FemInterface& FemL, const LocalVector& U, int ile, int col, LocalNodeData& Q, const LocalCellData& QC) const
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
      BasicIntegrator::universal_point(FemL,_UH,U);
      BasicIntegrator::universal_point(FemL,_QH,Q);
      FemL.x(x);
      FemL.normal(n);
      double  h = FemL.G();
      double  weight = IF->w(k)*h;
      BE.SetFemData(_QH);
      BE.pointboundary(h,_UH,x,n);
      for (int i=0;i<FemH.n();i++)
      {
        FemH.init_test_functions(_NN,weight,i);
        BE.Form(F.start(i),_UH,_NN,col);
      }
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void IntegratorQ1Q2<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, 
    const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalNodeData& Q, const LocalCellData& QC) const
{
  assert(FemH.n()==FemL.n());
  b.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x;
  FemH.point(p);
  FemL.point(p);
  FemL.x(x);
  BasicIntegrator::universal_point(FemL,_QH,Q);
  DRHS.SetFemData(_QH);

  for (int i=0; i<FemH.n(); i++)
    {
      FemH.init_test_functions(_NN,1.,i);
      DRHS.operator()(j,b.start(i),_NN,x);
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
double IntegratorQ1Q2<DIM>::IntegratorQ1Q2::MassMatrix(EntryMatrix& E, const FemInterface& FemH, 
    const FemInterface& FemL) const
{
  FemFunction NnnH, NnnL;

  NnnH.resize(FemH.n());
  NnnL.resize(FemL.n());
  
  E.SetDimensionDof(FemH.n(),FemL.n());
  int ncomp=1;
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();
  
  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;
  double omega = 0.;
  for (int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point(xi);
    FemL.point(xi);
    double vol = FemL.J();
    double weight  = IF->w(k) * vol;
    omega += weight;
    for (int i=0;i<FemH.n();i++)
    {
      FemH.init_test_functions(NnnH[i],1.,i);
    }
    for (int i=0;i<FemL.n();i++)
    {
      FemL.init_test_functions(NnnL[i],1.,i);
    }
    for (int i=0;i<FemH.n();i++)
    {
      for (int j=0;j<FemL.n();j++)
      {
        E.SetDofIndex(i,j);
        E(0,0) += weight * NnnH[i].m() * NnnL[j].m();
      }
    }
  }
  delete IF;
  return omega;
}

/*---------------------------------------------------------*/

template class IntegratorQ1Q2<2>;
template class IntegratorQ1Q2<3>;
}
