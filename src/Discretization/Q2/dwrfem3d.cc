#include "dwrfem.h" 
#include "galerkinintegratorq2.h"
#include "baseq23d.h"
#include "integratorq1q2.h"

namespace Gascoigne
{

/*---------------------------------------------------*/

DwrFem3d::DwrFem3d() : Q23d()
{
}

/*---------------------------------------------------*/

void DwrFem3d::BasicInit(const ParamFile*  paramfile)
{
  assert(PatchMeshInterpretor::GetIntegrator()==NULL);
  PatchMeshInterpretor::GetIntegratorPointer() =  new IntegratorQ1Q2<3>;

  assert(PatchMeshInterpretor::GetFem()==NULL);
  typedef Transformation3d<BaseQ23d>           TransQ2;
  typedef FiniteElement<3,2,TransQ2,BaseQ23d>  FiniteElement;

  PatchMeshInterpretor::GetFemPointer() =  new FiniteElement;
  PatchMeshInterpretor::BasicInit(paramfile);
}

/*---------------------------------------------------*/

void DwrFem3d::TransformationQ1(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  for(int ii=0;ii<ne;ii++)
    {
      Vertex3d v = GetMesh()->vertex3d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
      T(2,ii) = v.z();
    }
}

/*---------------------------------------------------*/

void DwrFem3d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->Form(EQ,__F,HighOrderFem,LowOrderFem,__U,__Q);
      PatchMeshInterpretor::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFem3d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(iq);
      I->Rhs(RHS,__F,HighOrderFem,LowOrderFem,__Q);
      PatchMeshInterpretor::LocalToGlobal(f,__F,iq,s);
    }
}
}
