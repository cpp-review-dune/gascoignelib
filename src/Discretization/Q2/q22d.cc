#include  "q22d.h"
#include  "galerkinintegratorq2.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq22d.h"
#include  "sparsestructure.h"
#include  "gascoignemeshtransfer.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "hnstructureq22d.h"

using namespace Gascoigne;
using namespace std;

/* ----------------------------------------- */

Q22d::Q22d() : Q2()
{
  HN = new HNStructureQ22d;
}

/* ----------------------------------------- */

Q22d::~Q22d()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

nmatrix<double> Q22d::GetLocalInterpolationWeights(int iq) const
{
  int nn = GetMesh()->nodes_per_cell(iq);
  nmatrix<double> w(nn,nn);
  w.zero();
  w(0,1) =  0.5  ; w(0,2) =  0.5  ; w(0,3) =  0.25;
  w(1,0) =  0.5  ; w(1,2) =  0.25 ; w(1,3) =  0.5 ;
  w(2,0) =  0.5  ; w(2,1) =  0.25 ; w(2,3) =  0.5 ;
  w(3,0) =  0.25 ; w(3,1) =  0.5  ; w(3,2) =  0.5 ;
  return w;
}

/* ----------------------------------------- */

void Q22d::BasicInit(const ParamFile* paramfile)
{
  if (GetIntegrator()==NULL)
    PatchMeshInterpretor::GetIntegratorPointer() =  new GalerkinIntegratorQ2<2>;

  assert(PatchMeshInterpretor::GetFem()==NULL);
  typedef Transformation2d<BaseQ22d>           TransQ2;
  typedef FiniteElement<2,1,TransQ2,BaseQ22d>  FiniteElement;

  PatchMeshInterpretor::GetFemPointer() =  new FiniteElement;
  PatchMeshInterpretor::BasicInit(paramfile);
}

/* ----------------------------------------- */

void Q22d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
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
  const GascoigneMeshTransfer* GMT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
  assert(GMT);
  //
  // dast ist einfach von Q12d kopiert !!!!!
  //

  const map<int,fixarray<2,int> >& zweier = GMT->GetZweier();
  const map<int,fixarray<4,int> >& vierer = GMT->GetVierer();
  const map<int,fixarray<8,int> >& achter = GMT->GetAchter();
  const nvector<int>& c2f    = GMT->GetC2f();

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
