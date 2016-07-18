/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#include  "q13d.h"
#include  "galerkinintegrator.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq1.h"
#include  "sparsestructure.h"
#include  "hnstructureq13d.h"
#include  "gascoignemesh.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"
#include  "hnstructureq12d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q13d::Q13d() : Q1() 
{
}

/* ----------------------------------------- */

/* ----------------------------------------- */

void Q13d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  if(!Q1<3>::GetIntegratorPointer())
    Q1<3>::GetIntegratorPointer() =  new GalerkinIntegrator<3>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  if(!Q1<3>::GetFemPointer())
    {
      typedef Transformation3d<BaseQ1<3> >           TransQ1;
      typedef FiniteElement<3,2,TransQ1,BaseQ1<3> >  FiniteElement;
      Q1<3>::GetFemPointer() =  new FiniteElement;
    }
  assert(GetFem());

  Q1<3>::BasicInit(pf);
}

/* ----------------------------------------- */

nmatrix<double> Q13d::GetLocalInterpolationWeights() const
{
  // w(i,j) = interpolation weight of node i to node j
  int nn = 8;//GetMesh()->nodes_per_cell();
  nmatrix<double> w(nn,nn);
  w.zero();

  w(0,1) =  0.5 ; w(0,2) = 0.5 ; w(0,3) = 0.25; w(0,4) = 0.5  ; w(0,5) = 0.25 ; w(0,6) = 0.25 ; w(0,7) = 0.125; 
  w(1,0) =  0.5 ; w(1,2) = 0.25; w(1,3) = 0.5 ; w(1,4) = 0.25 ; w(1,5) = 0.5  ; w(1,6) = 0.125; w(1,7) = 0.25;
  w(2,0) =  0.5 ; w(2,1) = 0.25; w(2,3) = 0.5 ; w(2,4) = 0.25 ; w(2,5) = 0.125; w(2,6) = 0.5  ; w(2,7) = 0.25;
  w(3,0) =  0.25; w(3,1) = 0.5 ; w(3,2) = 0.5 ; w(3,4) = 0.125; w(3,5) = 0.25 ; w(3,6) = 0.25 ; w(3,7) = 0.5;
  w(4,0) =  0.5 ; w(4,1) = 0.25; w(4,2) = 0.25; w(4,3) = 0.125; w(4,5) = 0.5  ; w(4,6) = 0.5  ; w(4,7) = 0.25; 
  w(5,0) =  0.25; w(5,1) = 0.5 ; w(5,2) = 0.125;w(5,3) = 0.25 ; w(5,4) = 0.5  ; w(5,6) = 0.25 ; w(5,7) = 0.5; 
  w(6,0) =  0.25; w(6,1) = 0.125;w(6,2) = 0.5 ; w(6,3) = 0.25 ; w(6,4) = 0.5  ; w(6,5) = 0.25 ; w(6,7) = 0.5; 
  w(7,0) =  0.125;w(7,1) = 0.25; w(7,2) = 0.25; w(7,3) = 0.5  ; w(7,4) = 0.25 ; w(7,5) = 0.5  ; w(7,6) = 0.5; 

  return w;
}


/* ----------------------------------------- */

void Q13d::StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const vector<int>& comp, double d) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  DoubleVector ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  FemData QH;

  for(int i=0;i<bv.size();i++)
  {
    int index = bv[i];
    const Vertex3d& v = GMP->vertex3d(index);

    QH.clear();
    GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
    for(; p!=GetDataContainer().GetNodeData().end(); p++)
    {
      QH[p->first].resize(p->second->ncomp());
      for(int c=0; c<p->second->ncomp(); c++)
      {
        QH[p->first][c].m() = p->second->operator()(index,c);
      }
    }

    BF.SetFemData(QH);

    BF(ff,v,col);
    for(int iii=0;iii<comp.size();iii++)
    {
      int c = comp[iii];
      u(index,c) = d * ff[c];
    }
  }
}

/* ----------------------------------------- */

void Q13d::Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
{
  if (&U==NULL) return;

  for(int in=0; in<GetMesh()->nnodes(); ++in)
    {
      Vertex3d v = GetMesh()->vertex3d(in);
      for(int c=0;c<u.ncomp();c++)
	{
	  u(in,c) = U(c,v);
	}
    }
}

/* ----------------------------------------- */

void Q13d::EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype, double d_visc) const
{
  EnergyEstimatorIntegrator<3> EEI(s_energytype,d_visc);
  const HierarchicalMesh3d*    HM = dynamic_cast<const HierarchicalMesh3d*>(EIC.GetMesh());

  EdgeInfoContainer<3>& EICC = dynamic_cast<EdgeInfoContainer<3>&>(EIC);

  EEI.BasicInit();

  // Kanten initialisieren
  EEJumps(EICC,u,EEI,HM);
  
  // Kantenintegrale auswerten
  EEJumpNorm(EICC,eta,EEI,HM);

  // Residuenterme auswerten
  EEResidual(eta,u,EQ,RHS,EEI);
}

/* ----------------------------------------- */

void Q13d::EEJumps(EdgeInfoContainer<3>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const
{
  fixarray<4,int> vertexes;
  nmatrix<double> T;

  for(int iq=0;iq<HM->ncells();++iq)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      GlobalToLocal_HM(__U,u,HM,iq);
      
      for (int ile=0; ile<6; ile++)
      {
        EEI.Jumps(__F,*GetFem(),__U,ile);

        int facenumber = HM->face_of_hex(iq,ile);

        if (EIC[facenumber]==NULL)
        {
          const Edge& edge = HM->edge(facenumber);
          HM->HexLawOrder().globalvertices_of_face(HM->hex(edge.master()),vertexes,edge.LocalMasterIndex());
          EIC[facenumber] = new EdgeInfo<3>();
          EIC[facenumber]->BasicInit(&edge,u.ncomp(),vertexes);
        }
        EIC[facenumber]->AddNodes(__F);
      }
    }
  }

  EIC.ModifyHanging();
}

/* ----------------------------------------- */

void Q13d::EEJumpNorm(EdgeInfoContainer<3>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<3>& EEI, const HierarchicalMesh3d* HM) const
{
  nmatrix<double> T;

  for (int iq=0; iq<HM->ncells(); iq++)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      double jump = 0.;
      for (int ile=0; ile<6; ile++)
      {
        int facenumber = HM->face_of_hex(iq,ile);
        if (EIC[facenumber]->GetCount()==2)
        {
          jump += EEI.JumpNorm(*GetFem(),EIC[facenumber]->GetNorm(),ile);
        }
      }
      for (int in=0; in<8; in++)
      {
        eta[HM->vertex_of_cell(iq,in)] += 0.125 * 0.5 * sqrt(jump);
      }
    }
  }
}

/* ----------------------------------------- */

void Q13d::EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const EnergyEstimatorIntegrator<3>& EEI) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  if (RHS) RHS->SetParameterData(__QP);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);

    //    GlobalToLocalData(iq);
    GlobalToLocal(__U,u,iq);
//    EQ.SetCellData(__QC);
    double res = EEI.Residual(__U,*GetFem(),EQ,RHS,__QN);
    for (int in=0; in<8; in++)
    {
      eta[GetMesh()->vertex_of_cell(iq,in)] += 0.125 * sqrt(res);
    }
  }
}

/* ----------------------------------------- */

int Q13d::GetCellNumber(const Vertex3d& p0, Vertex3d& p, int c0) const
{
  int iq;
  
  for(iq=0; iq<GetMesh()->ncells(); ++iq)
  {
    bool found = true;

    for(int d=0; d<3; ++d)
    {
      double min=GetMesh()->vertex3d(GetMesh()->vertex_of_cell(iq,0))[d];
      double max=min;
      for(int j=1; j<8; ++j)
      {
        double x = GetMesh()->vertex3d(GetMesh()->vertex_of_cell(iq,j))[d];
        
        min = Gascoigne::min(min,x);
        max = Gascoigne::max(max,x);
      }
      if((p0[d]<min)||(p0[d]>max)) 
      {
        found = false;
        break;
      }
    }
    
    if(!found)
    {
      continue;
    }
    
    VertexTransformation(p0,p,iq);
    
    for(int d=0; d<3; ++d)
    {
      if((p[d]<0.-1.e-12)||(p[d]>1.+1.e-12))
      {
        found = false;
      }
    }
    
    if(found)
    {
      break;
    }
  }

  if(iq<GetMesh()->ncells())
  {
    return iq;
  }
  else
  {
    return -1;
  }
}

/* ----------------------------------------- */

void Q13d::VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation3d<BaseQ1<3> > Tr;
  Tr.init(T);

  Vertex3d res;
  
  p = 0.5;
  
  for(int niter=1; ;niter++)
  {
    Tr.point(p);
    
    res = p0;
    res.add(-1,Tr.x());

    if(res.norm()<1.e-13)
    {
      break;
    }
    assert(niter<10);
    
    Tr.DTI().mult_ad(p,res);
  } 
}

/* ----------------------------------------- */

}
