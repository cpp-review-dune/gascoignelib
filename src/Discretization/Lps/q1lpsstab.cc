#include  "q1lpsstab.h"
#include  "transformation2d.h"
#include  "transformation3d.h"
#include  "finiteelement.h"
#include  "baseq1patch.h"
#include  "baseq13dpatch.h"
#include  "lpsintegrator.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q1LpsStab::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  PatchMeshInterpretor::BasicInit(paramfile);
  HN = hn;
  assert(HN);
}

/* ----------------------------------------- */

void Q1LpsStab::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  HN->CondenseHangingPatch(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}

/* ----------------------------------------- */

void Q1LpsStab2d::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  assert(PatchMeshInterpretor::GetIntegrator()==NULL);
  PatchMeshInterpretor::GetIntegratorPointer() =  new LpsIntegratorQ1<2>;
  
  assert(PatchMeshInterpretor::GetFem()==NULL);
  typedef Transformation2d<BaseQ12d>           TransQ1;
  typedef FiniteElement<2,1,TransQ1,BaseQ12dPatch>  FiniteElement;
  PatchMeshInterpretor::GetFemPointer() =  new FiniteElement;
  
  Q1LpsStab::BasicInit(paramfile,hn);  
}

/* ----------------------------------------- */

void Q1LpsStab3d::BasicInit(const ParamFile* paramfile, const HNStructureInterface* hn)
{
  assert(PatchMeshInterpretor::GetIntegrator()==NULL);
  PatchMeshInterpretor::GetIntegratorPointer() =  new LpsIntegratorQ1<3>;
  
  assert(PatchMeshInterpretor::GetFem()==NULL);
  typedef Transformation3d<BaseQ13d>           TransQ1;
  typedef FiniteElement<3,2,TransQ1,BaseQ13dPatch>  FiniteElement;
  PatchMeshInterpretor::GetFemPointer() =  new FiniteElement;
  
  Q1LpsStab::BasicInit(paramfile,hn);  
}

/* ----------------------------------------- */

void Q1LpsStab2d::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = 2;
  int ne = GetPatchMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);

  for(int ii=0;ii<ne;ii++)
    {
      Vertex2d v = GetPatchMesh()->vertex2d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
    }
}

/* ----------------------------------------- */

void Q1LpsStab3d::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = 3;
  int ne = GetPatchMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);

  for(int ii=0;ii<ne;ii++)
    {
      Vertex3d v = GetPatchMesh()->vertex3d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
      T(2,ii) = v.z();
    }
}

}