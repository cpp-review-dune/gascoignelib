#include "dwr_integrator.h"
#include "patchintegrationformula.h"

namespace Gascoigne
{
/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::BasicInit()
{
}


/*---------------------------------------------------*/

/**
 * Wandelt die lokale Nummerierung der Zellen eines Patches von der Gascoigne-nummerierung in die 
 * vom Integrator verwendete um! Funktioniert aus Symmetrie-Gruenden auch andersherum.
 * Erwartet eine Zahl zwischen 0 und 3 (2-D) bzw 0 bis 7 (3-D)
*/
  
  template<int DIM>
  int DWR_Integrator<DIM>::PatchMeshNr2IntegratorNr(int in) const
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
  
  template<int DIM>
  int DWR_Integrator<DIM>::PatchMeshNr2IntegratorNrBoundary(int in, int ile) const
  { 
    if(DIM==2)
      {
	if(ile == 0) return in;
	if(ile == 1) return in+1;
	if(ile == 2) return 3-in;
	if(ile == 3) return 3*in;
	std::cerr<<"Integrator Q1/Q2: Invalid ile "<<ile<<"!"<<std::endl;
	abort();
      }
    if(DIM==3)
      {
	if(ile == 0)
	  {
	    if(in == 3) return 2;
	    if(in == 2) return 3;
	    return in;
	  } 
	if(ile == 1)
	  {
	    if(in == 0) return 1;
	    if(in == 1) return 5;
	    if(in == 2) return 2;
	    return 6;
	  }
	if(ile == 2)
	  {
	    if(in == 0) return 3;
	    if(in == 1) return 2;
	    if(in == 2) return 7;
	    return 6;
	  }
	if(ile == 3)
	  {
	    if(in == 0) return 4;
	    if(in == 1) return 0;
	    if(in == 2) return 7;
	    return 3;
	  }
	if(ile == 4)
	  {
	    if(in == 0) return 4;
	    if(in == 1) return 5;
	    if(in == 2) return 0;
	    return 1;
	  }
	if(ile == 5)
	  {
	    if(in == 0) return 5;
	    if(in == 1) return 4;
	    if(in == 2) return 6;
	    return 7;
	  }
      }
    std::cerr<<"Integrator Q1/Q2: Dimension "<<DIM<<" nicht implementiert!"<<std::endl;
    abort();
  }

/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::Form(const Equation& EQ, LocalVector& F, 
			       const FemInterface& FemAnsatz, const FemInterface& FemTrial, 
			       const LocalVector& U, const LocalData& Q, const LocalData& QC) const
{
  //  assert(FemAnsatz.n()==FemTrial.n());

  F.ReInit(U.ncomp(),FemTrial.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;

  int numcells;
  if(DIM==2) numcells = 4;
  else       numcells = 8;
  
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
      FemTrial.point(xi);
      FemAnsatz.point(xi);
      double vol = FemTrial.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemAnsatz,_UH,U);
      BasicIntegrator::universal_point(FemAnsatz,_QH,Q);
      
      FemAnsatz.x(x);

      EQ.SetFemData(_QH);
      //Die Zelldatensetzen
      if(k%numcells == 0)
      {
	  //Wir sind am anfang einer neuen Zelle
	  //Die CellData aus den LocalData hohlen.
	  int IntegrCellNum = k/numcells;
	  int PatchCellNum=PatchMeshNr2IntegratorNr(IntegrCellNum);
	  
	  universal_point(_QCH,QC,PatchCellNum);
	  EQ.SetCellData(_QCH);
      }
      
      EQ.point(h,_UH,x);

      for (int i=0;i<FemTrial.n();i++)
	{
	  FemTrial.init_test_functions(_NN,weight,i);
	  EQ.Form(F.start(i),_UH,_NN);
	}
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, 
				      const FemInterface& FemAnsatz, const FemInterface& FemTrial, 
				      const LocalVector& Z, const LocalData& Q, const LocalData& QC) const
{
  assert(FemAnsatz.n()==FemTrial.n());

  F.ReInit(Z.ncomp(),FemTrial.n());
  F.zero();

  _NNN.resize(FemTrial.n());
  EntryMatrix E;
  E.SetDimensionDof(FemTrial.n(),FemAnsatz.n());
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
      FemAnsatz.point(xi);
      FemTrial.point(xi);
      double vol = FemTrial.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemAnsatz,_UH,Z);
      BasicIntegrator::universal_point(FemAnsatz,_QH,Q);
      FemAnsatz.x(x);
      //EQ.pointmatrix(h,_QH["u"],_QH,x);
      EQ.pointmatrix(h,_QH["u"],x);
      for (int i=0;i<FemTrial.n();i++)
	{
	  FemTrial.init_test_functions(_NNN[i],weight,i);
	}
      for (int j=0;j<FemAnsatz.n();j++)
	{
	  FemAnsatz.init_test_functions(MM,1.,j);
	  for (int i=0;i<FemTrial.n();i++)
	    {
	      E.SetDofIndex(j,i);
	      EQ.Matrix(E,_QH["u"],_NNN[i],MM);
	    }
	}
    }
  for (int i=0;i<FemTrial.n();i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  double sum = 0.;
	  for (int j=0; j<FemAnsatz.n(); j++)
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
void DWR_Integrator<DIM>::Rhs(const DomainRightHandSide& RHS, LocalVector& F, 
			      const FemInterface& FemAnsatz, const FemInterface& FemTrial, 
			      const LocalData& Q, const LocalData& QC) const
{
  //assert(FemAnsatz.n()==FemTrial.n());
  
  F.ReInit(RHS.GetNcomp(),FemTrial.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
  else        IF = new PatchFormula3d<27,HexGauss27>;

  Vertex<DIM> x, xi;

  int numcells;
  if(DIM==2) numcells = 4;
  else       numcells = 8;
  
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
      FemAnsatz.point(xi);
      FemTrial.point(xi);
      double vol = FemAnsatz.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemAnsatz,_QH,Q);
      RHS.SetFemData(_QH);
      //Die Zelldatensetzen
      if(k%numcells == 0)
	{
	  int IntegrCellNum = k/numcells;
	  int PatchCellNum=PatchMeshNr2IntegratorNr(IntegrCellNum);
	  
	  universal_point(_QCH,QC,PatchCellNum);
	  RHS.SetCellData(_QCH);
      }
      
      RHS.SetCellSize(h);
      FemTrial.x(x);

      for (int i=0;i<FemTrial.n();i++)
	{
	  FemTrial.init_test_functions(_NN,weight,i);
	  RHS(F.start(i),_NN,x);
	}
    } 
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, 
				      const FemInterface& FemAnsatz, const FemInterface& FemTrial,
				      int ile, int col, const LocalData& Q, const LocalData& QC) const
{
  //  hier nochmal ueberlegen. wie muss der rand approximiert werden, z.b. welcher normalvektor?
  abort();
  
//   assert(FemAnsatz.n()==FemTrial.n());

//   F.ReInit(RHS.GetNcomp(),FemAnsatz.n());
//   F.zero();

//   IntegrationFormulaInterface* IF;
//   if (DIM==2) IF = new PatchFormula1d<3,LineGauss3>;
//   else        IF = new PatchFormula2d<9,QuadGauss9>;

//   Vertex<DIM> x, n;
//   Vertex<DIM-1> xi;

//   int numcells;
//   if(DIM==2) numcells = 2;
//   else numcells = 4;
  
//   //Soviele Integr. Pkte je Zelle
//   numcells = IF->n()/numcells;
//   //Test, ist das sinnvollgewesen
//   if(IF->n()%numcells != 0)
//   {
//       std::cerr<<"Integrator Q1/Q2: Falsche Anzahl Zellen je Patch oder unzulaessige Integrationsformel!"<<std::endl;
//       abort();
//   }

//   for (int k=0; k<IF->n(); k++)
//     {
//       IF->xi(xi,k);
//       FemAnsatz.point_boundary(ile,xi);
//       FemTrial.point_boundary(ile,xi);
//       BasicIntegrator::universal_point(FemL,_QH,Q);
//       RHS.SetFemData(_QH);
//       FemTrial.x(x);
//       FemTrial.normal(n);
//       double  h = FemL.G();
//       double  weight = IF->w(k)*h;
//       //Die Zelldatensetzen
//       if(k%numcells == 0)
//       {
// 	  //Wir sind am anfang einer neuen Zelle
// 	  //Die CellData aus den LocalData hohlen.
// 	  int IntegrCellNum = k/numcells;
// 	  int PatchCellNum=PatchMeshNr2IntegratorNrBoundary(IntegrCellNum,ile);
	  
// 	  universal_point(_QCH,QC,PatchCellNum);
// 	  RHS.SetCellData(_QCH);
//       }      
//       RHS.SetCellSize(h);
//       for (int i=0;i<FemH.n();i++)
//       {
//         FemH.init_test_functions(_NN,weight,i);
//         RHS(F.start(i),_NN,x,n,col);
//       }
//     }
//   delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemAnsatz, 
    const FemInterface& FemTrial, const LocalVector& U, int ile, int col, LocalData& Q, const LocalData& QC) const
{
  F.ReInit(U.ncomp(), FemTrial.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula1d<3,LineGauss3>;
  else        IF = new PatchFormula2d<9,QuadGauss9>;
  
  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;
  
  int numcells;
  if(DIM==2) numcells = 2;
  else numcells = 4;
  
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
      FemTrial.point_boundary(ile,xi);
      FemAnsatz.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FemAnsatz,_UH,U);
      BasicIntegrator::universal_point(FemAnsatz,_QH,Q);

      FemAnsatz.x(x);
      FemAnsatz.normal(n);
      
      double  h =  FemTrial.G();
      double  weight = IF->w(k)*h;
      BE.SetFemData(_QH);

      //Die Zelldatensetzen
      int IntegrCellNum = -1;
      if(k%numcells == 0)
      {
	  //Wir sind am anfang einer neuen Zelle
	  //Die CellData aus den LocalData hohlen.
	  ///int IntegrCellNum = k/numcells;
	  IntegrCellNum = k/numcells;
	  int PatchCellNum=PatchMeshNr2IntegratorNrBoundary(IntegrCellNum,ile);
	 
	  BasicIntegrator::universal_point(_QCH,QC,PatchCellNum);
	  BE.SetCellData(_QCH);
      }
      BE.pointboundary(h,_UH,x,n);
      for (int i=0;i<FemTrial.n();i++)
      {
        FemTrial.init_test_functions(_NN,weight,i);
        BE.Form(F.start(i),_UH,_NN,col);
      }
    }
  delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
void DWR_Integrator<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, 
    const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalData& Q, const LocalData& QC) const
{
  abort();
  
//   assert(FemH.n()==FemL.n());
//   b.zero();

//   IntegrationFormulaInterface* IF;
//   if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
//   else        IF = new PatchFormula3d<27,HexGauss27>;

//   Vertex<DIM> x;
//   FemH.point(p);
//   FemL.point(p);
//   FemL.x(x);
//   BasicIntegrator::universal_point(FemL,_QH,Q);
//   DRHS.SetFemData(_QH);

//   for (int i=0; i<FemH.n(); i++)
//     {
//       FemH.init_test_functions(_NN,1.,i);
//       DRHS.operator()(j,b.start(i),_NN,x);
//     }
//   delete IF;
}

/*---------------------------------------------------*/

template<int DIM>
double DWR_Integrator<DIM>::DWR_Integrator::MassMatrix(EntryMatrix& E, const FemInterface& FemH, 
    const FemInterface& FemL) const
{
  abort();
  
//   FemFunction NnnH, NnnL;

//   NnnH.resize(FemH.n());
//   NnnL.resize(FemL.n());
  
//   E.SetDimensionDof(FemH.n(),FemL.n());
//   int ncomp=1;
//   E.SetDimensionComp(ncomp,ncomp);
//   E.resize();
//   E.zero();
  
//   IntegrationFormulaInterface* IF;
//   if (DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
//   else        IF = new PatchFormula3d<27,HexGauss27>;

//   Vertex<DIM> x, xi;
//   double omega = 0.;
//   for (int k=0; k<IF->n(); k++)
//   {
//     IF->xi(xi,k);
//     FemH.point(xi);
//     FemL.point(xi);
//     double vol = FemL.J();
//     double weight  = IF->w(k) * vol;
//     omega += weight;
//     for (int i=0;i<FemH.n();i++)
//     {
//       FemH.init_test_functions(NnnH[i],1.,i);
//     }
//     for (int i=0;i<FemL.n();i++)
//     {
//       FemL.init_test_functions(NnnL[i],1.,i);
//     }
//     for (int i=0;i<FemH.n();i++)
//     {
//       for (int j=0;j<FemL.n();j++)
//       {
//         E.SetDofIndex(i,j);
//         E(0,0) += weight * NnnH[i].m() * NnnL[j].m();
//       }
//     }
//   }
//   delete IF;
//   return omega;
}

/*---------------------------------------------------------*/

template<int DIM>
void Gascoigne::DWR_Integrator<DIM>::MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U) const
{
  abort();
  
//   assert(FemH.n()==FemL.n());

//   F.ReInit(U.ncomp(),FemH.n());
//   F.zero();

//   IntegrationFormulaInterface* IF;
//   if(DIM==2) IF = new PatchFormula2d<9,QuadGauss9>;
//   else       IF = new PatchFormula3d<27,HexGauss27>;

//   Vertex<DIM> xi;

//   for(int k=0; k<IF->n(); k++)
//   {
//     IF->xi(xi,k);
//     FemH.point(xi);
//     FemL.point(xi);
//     double vol = FemL.J();
//     double weight = IF->w(k) * vol;
//     BasicIntegrator::universal_point(FemL,_UH,U);
//     for(int i=0;i<FemH.n();i++)
//     {
//       FemH.init_test_functions(_NN,weight,i);
//       for(int m=0; m<TP.n(); m++)
//       {
//         for(int n=0; n<TP.n(); n++)
//         {
//           F(i,m) += TP(m,n) * _UH[n].m() * _NN.m();
//         }
//       }
//     }
//   }
//   delete IF;
}

/*---------------------------------------------------------*/

template class DWR_Integrator<2>;
template class DWR_Integrator<3>;
}
