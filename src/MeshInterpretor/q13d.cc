#include  "q13d.h"
#include  "galerkinintegrator.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq13d.h"
#include  "sparsestructure.h"
#include  "hnstructureq13d.h"
#include  "gascoignemesh.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"

using namespace std;

/* ----------------------------------------- */

Q13d::Q13d() : Q1() 
{
  HN = new HNStructureQ13d;
}

/* ----------------------------------------- */

Q13d::~Q13d()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

void Q13d::BasicInit(const std::string& paramfile)
{
  assert(CellMeshInterpretor::GetIntegrator()==NULL);
  BasicMeshInterpretor::GetIntegratorPointer() =  new GalerkinIntegrator<3>;

  assert(CellMeshInterpretor::GetFem()==NULL);
  typedef Transformation3d<BaseQ13d>           TransQ1;
  typedef FiniteElement<3,2,TransQ1,BaseQ13d>  FiniteElement;
  BasicMeshInterpretor::GetFemPointer() =  new FiniteElement;

  CellMeshInterpretor::BasicInit(paramfile);
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

void Q13d::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  nvector<double> ff(u.ncomp(),0.);
  const IntVector& bv = GMP->VertexOnBoundary(col);

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];
      const Vertex3d& v = GMP->vertex3d(index);
      
      BF(ff,v,col);
      for(int iii=0;iii<comp.size();iii++)
	{
	  int c = comp[iii];
	  u(index,c) = ff[c];
	}
    }
}

/* ----------------------------------------- */

void Q13d::Interpolate(GlobalVector& u, const InitialCondition& U) const
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

void Q13d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
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
      SS.build_add(i,c2f[i]);
    }
  for(std::map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      for(int ii=0;ii<2;ii++) SS.build_add(il,n2[ii]);
    }
  for(std::map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    for(int ii=0;ii<4;ii++) SS.build_add(il,n4[ii]);
  }
  for(std::map<int,fixarray<8,int> >::const_iterator p=achter.begin();
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
      val[ST.Find(c2f[i],i)] = 1.;
    }
  for(std::map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      val[ST.Find(il,n2[0])] = 0.5;
      val[ST.Find(il,n2[1])] = 0.5;
    }
  for(std::map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    val[ST.Find(il,n4[0])] = 0.25;
    val[ST.Find(il,n4[1])] = 0.25;
    val[ST.Find(il,n4[2])] = 0.25;
    val[ST.Find(il,n4[3])] = 0.25;
  }
  for(std::map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	val[ST.Find(il,n8[i])] = 0.125;
      }
  }
}
