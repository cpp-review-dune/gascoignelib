#include  "q23d.h"
#include  "galerkinintegratorq2.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq23d.h"
#include  "sparsestructure.h"
#include  "gascoignemeshtransfer.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "hnstructureq23d.h"

using namespace Gascoigne;
using namespace std;

/* ----------------------------------------- */

Q23d::Q23d() : Q2()
{
  HN = new HNStructureQ23d;
}

/* ----------------------------------------- */

Q23d::~Q23d()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

nmatrix<double> Q23d::GetLocalInterpolationWeights(int iq) const
{
  int nn = GetMesh()->nodes_per_cell(iq);
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

void Q23d::BasicInit(const ParamFile* paramfile)
{
  assert(PatchMeshInterpretor::GetIntegrator()==NULL);
  PatchMeshInterpretor::GetIntegratorPointer() =  new GalerkinIntegratorQ2<3>;

  assert(PatchMeshInterpretor::GetFem()==NULL);
  typedef Transformation3d<BaseQ23d>           TransQ2;
  typedef FiniteElement<3,2,TransQ2,BaseQ23d>  FiniteElement;

  PatchMeshInterpretor::GetFemPointer() =  new FiniteElement;
  PatchMeshInterpretor::BasicInit(paramfile);
}

/* ----------------------------------------- */

void Q23d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    if(IP)
      {
	IP->BasicInit(MT);
	return;
      }
  }
  assert(0);
}
