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
using namespace Gascoigne;

/* ----------------------------------------- */

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

  assert(CellMeshInterpretor::GetIntegrator()==NULL);
  BasicMeshInterpretor::GetIntegratorPointer() =  new GalerkinIntegrator<2>;

  assert(CellMeshInterpretor::GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
  BasicMeshInterpretor::GetFemPointer() =  new FiniteElement;

  CellMeshInterpretor::BasicInit(pf);
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
  nvector<double> ff(u.ncomp(),0.);
  const IntVector& bv = GMP->VertexOnBoundary(col);

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

void Q12d::Interpolate(GlobalVector& u, const InitialCondition& U) const
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

void Q12d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    if(IP)
      {
	IP->init(MT);
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
  const nvector<int>& c2f    = GT->GetC2f();

  int n  = c2f.size() +   zweier.size() +   vierer.size() +   achter.size();
  int nt = c2f.size() + 2*zweier.size() + 4*vierer.size() + 8*achter.size();

  ColumnStencil& ST = IP->GetStencil();
  nvector<double>& val = IP->GetAlpha();

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
