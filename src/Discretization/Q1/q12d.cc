#include  "q12d.h"
#include  "galerkinintegrator.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"
#include  "sparsestructure.h"
#include  "gascoignemesh.h"
#include  "hnstructureq12d.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"
#include  "hnstructureq12d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q12d::Q12d() : Q1()
{
}

/* ----------------------------------------- */

HNStructureInterface* Q12d::NewHNStructure()
{
  return new HNStructureQ12d;
}

/* ----------------------------------------- */

void Q12d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  if(!GetIntegratorPointer())
    GetIntegratorPointer() =  new GalerkinIntegrator<2>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  if(!GetFemPointer())
    {
      typedef Transformation2d<BaseQ12d>           TransQ1;
      typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
      CellDiscretization::GetFemPointer() =  new FiniteElement;
    }
  assert(GetFem());

  CellDiscretization::BasicInit(pf);
}

/* ----------------------------------------- */

nmatrix<double> Q12d::GetLocalInterpolationWeights() const
{
  // w(i,j) = interpolation weight of node i to node j
  int nn = 4;//GetMesh()->nodes_per_cell();
  nmatrix<double> w(nn,nn);
  w.zero();
  w(0,1) =  0.5  ; w(0,2) =  0.5  ; w(0,3) =  0.25;
  w(1,0) =  0.5  ; w(1,2) =  0.25 ; w(1,3) =  0.5 ;
  w(2,0) =  0.5  ; w(2,1) =  0.25 ; w(2,3) =  0.5 ;
  w(3,0) =  0.25 ; w(3,1) =  0.5  ; w(3,2) =  0.5 ;
  return w;
}

/* ----------------------------------------- */

void Q12d::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  DoubleVector ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  GlobalToGlobalData();
  BF.SetParameterData(__qq);

  for(int ii=0;ii<comp.size();ii++)
    {
      int c = comp[ii];
      if(c<0) {
        cerr << "negative component: " << c << endl;
        assert(0);
      } else if(c>=u.ncomp()){
        cerr << "unknown component: " << c << endl;
        assert(0);
      }
    }

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];
      const Vertex2d& v = GMP->vertex2d(index);
      
      BF(ff,v,col);
      for(int iii=0;iii<comp.size();iii++)
        {
          int c = comp[iii];
          u(index,c) = ff[c];
        }
    }
}

/* ----------------------------------------- */

void Q12d::Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
{
  if (&U==NULL) return;

  for(int in=0; in<GetMesh()->nnodes(); ++in)
    {
      Vertex2d v = GetMesh()->vertex2d(in);
      for(int c=0;c<u.ncomp();c++)
        {
          u(in,c) = U(c,v);
        }
    }
}

/* ----------------------------------------- */

void Q12d::InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const
{
  const IntVector& vo2n = *GetMesh()->Vertexo2n();
  nvector<bool> habschon(GetMesh()->nnodes(),0);  

  assert(vo2n.size()==uold.n());
  assert(GetMesh()->nnodes()==u.n());
  assert(u.ncomp()==uold.ncomp());

  for(int i=0;i<vo2n.size();i++)
    {
      int in = vo2n[i];

      if(in>=0) 
        {
          u.equ_node(in,1.,i,uold);
          habschon[in] = 1;
        }
    }
  nvector<fixarray<3,int> > nodes(4);
  nodes[0][0] = 1; nodes[0][1] = 0;  nodes[0][2] = 2;
  nodes[1][0] = 3; nodes[1][1] = 0;  nodes[1][2] = 6;
  nodes[2][0] = 5; nodes[2][1] = 2;  nodes[2][2] = 8;
  nodes[3][0] = 7; nodes[3][1] = 6;  nodes[3][2] = 8;
 
  const PatchMesh* PM = dynamic_cast<const PatchMesh*>(GetMesh());
  assert(PM);

  for(int iq=0;iq<PM->npatches();++iq)
    {
      IntVector vi =* PM->IndicesOfPatch(iq);

      for(int j=0; j<nodes.size(); j++)
        {
          int v  = vi[nodes[j][0]];
          int v1 = vi[nodes[j][1]];
          int v2 = vi[nodes[j][2]];
          assert(habschon[v1]);
          assert(habschon[v2]);
          if (habschon[v]==0) 
            {
              u.equ_node(v,0.5,v1,uold);
              u.add_node(v,0.5,v2,uold);
              habschon[v] = 1;
            }
        }
      int v = vi[4];
      if (habschon[v]==0)
        {
          u.equ_node(v,0.25,vi[0],uold);
          u.add_node(v,0.25,vi[2],uold);	  
          u.add_node(v,0.25,vi[6],uold);	  
          u.add_node(v,0.25,vi[8],uold);	  
          habschon[v] = 1;
        }
    }  
}

/* ----------------------------------------- */

void Q12d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    if(IP)
      {
        IP->BasicInit(MT);
        return;
      }
  }

  MgInterpolatorMatrix* IP = dynamic_cast<MgInterpolatorMatrix*>(I);
  assert(IP);
  const GascoigneMeshTransfer* GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
  assert(GT);

  const map<int,fixarray<2,int> >& zweier = GT->GetZweier();
  const map<int,fixarray<4,int> >& vierer = GT->GetVierer();
  const map<int,fixarray<8,int> >& achter = GT->GetAchter();
  const IntVector& c2f    = GT->GetC2f();

  int n  = c2f.size() +   zweier.size() +   vierer.size() +   achter.size();
  int nt = c2f.size() + 2*zweier.size() + 4*vierer.size() + 8*achter.size();

  ColumnStencil& ST = IP->GetStencil();
  DoubleVector& val = IP->GetAlpha();

  SparseStructure SS;

  SS.build_begin(n);
  for(int i=0;i<c2f.size();i++)
    {
      assert(c2f[i]>=0);

      SS.build_add(c2f[i],i);
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      for(int ii=0;ii<2;ii++) SS.build_add(il,n2[ii]);
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    for(int ii=0;ii<4;ii++) SS.build_add(il,n4[ii]);
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for(int ii=0;ii<8;ii++) SS.build_add(il,n8[ii]);
  }
  SS.build_end();

  assert(nt==SS.ntotal());

  ST.memory(&SS);

  val.reservesize(nt);

  for(int i=0;i<c2f.size();i++)
    {
      // ich weiss nicht, ob das richtig ist !!!!!
      int pos = ST.Find(c2f[i],i);
      assert(pos>=0);
      
      val[pos] = 1.;
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      val[ST.Find(il,n2[0])] = 0.5;
      val[ST.Find(il,n2[1])] = 0.5;
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    val[ST.Find(il,n4[0])] = 0.25;
    val[ST.Find(il,n4[1])] = 0.25;
    val[ST.Find(il,n4[2])] = 0.25;
    val[ST.Find(il,n4[3])] = 0.25;
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	val[ST.Find(il,n8[i])] = 0.125;
      }
  }
}

/* ----------------------------------------- */

void Q12d::EnergyEstimator(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS) const
{
  EnergyEstimatorIntegrator<2> EEI;
  const HierarchicalMesh2d*    HM = dynamic_cast<const HierarchicalMesh2d*>(EIC.GetMesh());

  EEI.BasicInit();

  // Kanten initialisieren
  EEJumps(EIC,u,EEI,HM);
  
  // Kantenintegrale auswerten
  EEJumpNorm(EIC,eta,EEI,HM);

  // Residuenterme auswerten
  EEResidual(eta,u,EQ,RHS,EEI);
}

/* ----------------------------------------- */

void Q12d::EnergyEstimatorZeroRhs(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ) const
{
  EnergyEstimatorIntegrator<2> EEI;
  const HierarchicalMesh2d*    HM = dynamic_cast<const HierarchicalMesh2d*>(EIC.GetMesh());

  EEI.BasicInit();

  // Kanten initialisieren
  EEJumps(EIC,u,EEI,HM);
  
  // Kantenintegrale auswerten
  EEJumpNorm(EIC,eta,EEI,HM);

  // Residuenterme auswerten
  EEResidualZeroRhs(eta,u,EQ,EEI);
}

/* ----------------------------------------- */

void Q12d::EEJumps(EdgeInfoContainer<2>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const
{
  fixarray<2,int> vertexes;
  nmatrix<double> T;

  for(int iq=0;iq<HM->ncells();++iq)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      GlobalToLocal_HM(__U,u,HM,iq);
      
      for (int ile=0; ile<4; ile++)
      {
        EEI.Jumps(__F,*GetFem(),__U,ile);

        int edgenumber = HM->edge_of_quad(iq,ile);

        if (EIC[edgenumber]==NULL)
        {
          const Edge& edge = HM->edge(edgenumber);
          HM->QuadLawOrder().globalvertices_of_edge(HM->quad(edge.master()),vertexes,edge.LocalMasterIndex());
          EIC[edgenumber] = new EdgeInfo<2>();
          EIC[edgenumber]->BasicInit(&edge,u.ncomp(),vertexes);
        }
        EIC[edgenumber]->AddNodes(__F);
      }
    }
  }

  EIC.ModifyHanging();
}

/* ----------------------------------------- */

void Q12d::EEJumpNorm(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const
{
  nmatrix<double> T;

  for (int iq=0; iq<HM->ncells(); iq++)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      double jump = 0.;
      for (int ile=0; ile<4; ile++)
      {
        int edgenumber = HM->edge_of_quad(iq,ile);
        if (EIC[edgenumber]->GetCount()==2)
        {
          jump += EEI.JumpNorm(*GetFem(),EIC[edgenumber]->GetNorm(),ile);
        }
      }
      for (int in=0; in<4; in++)
      {
        eta[HM->vertex_of_cell(iq,in)] += 0.25 * 0.5 * sqrt(jump);
      }
    }
  }
}

/* ----------------------------------------- */

void Q12d::EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide& RHS, const EnergyEstimatorIntegrator<2>& EEI) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  RHS.SetParameterData(__qq);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);
        
    GlobalToLocalData(iq);
    GlobalToLocal(__U,u,iq);
    double res = EEI.Residual(__U,*GetFem(),EQ,RHS,__Q);
    for (int in=0; in<4; in++)
    {
      eta[GetMesh()->vertex_of_cell(iq,in)] += 0.25 * sqrt(res);
    }
  }
}

/* ----------------------------------------- */

void Q12d::EEResidualZeroRhs(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const EnergyEstimatorIntegrator<2>& EEI) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);
        
    GlobalToLocalData(iq);
    GlobalToLocal(__U,u,iq);
    double res = EEI.ResidualZeroRhs(__U,*GetFem(),EQ,__Q);
    for (int in=0; in<4; in++)
    {
      eta[GetMesh()->vertex_of_cell(iq,in)] += 0.25 * sqrt(res);
    }
  }
}

/* ----------------------------------------- */

int Q12d::GetCellNumber(const Vertex2d& p0, Vertex2d& p) const
{
  int iq;
  
  for(iq=0; iq<GetMesh()->ncells(); ++iq)
  {
    bool found = true;

    for(int d=0; d<2; ++d)
    {
      double min=GetMesh()->vertex2d(GetMesh()->vertex_of_cell(iq,0))[d];
      double max=min;
      for(int j=1; j<4; ++j)
      {
        double x = GetMesh()->vertex2d(GetMesh()->vertex_of_cell(iq,j))[d];
        
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
    
    for(int d=0; d<2; ++d)
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

void Q12d::VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation2d<BaseQ12d> Tr;
  Tr.init(T);

  Vertex2d res;
  
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
