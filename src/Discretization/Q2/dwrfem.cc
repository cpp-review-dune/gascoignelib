#include "dwrfem.h" 
#include "galerkinintegratorq2.h"
#include "baseq22d.h"
#include "integratorq1q2.h"

namespace Gascoigne
{

/*---------------------------------------------------*/

DwrFem2d::DwrFem2d() : Q22d() {}

/*---------------------------------------------------*/

void DwrFem2d::BasicInit(const ParamFile*  paramfile)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() =  new IntegratorQ1Q2<2>;

  GetIntegratorPointer()->BasicInit();

  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation2d<BaseQ22d>           TransQ2;
  typedef FiniteElement<2,1,TransQ2,BaseQ22d>  FiniteElement;

  PatchDiscretization::GetFemPointer() =  new FiniteElement;
  PatchDiscretization::BasicInit(paramfile);
}

/*---------------------------------------------------*/

void DwrFem2d::TransformationQ1(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  nvector<int> indices = GetPatchMesh()->CoarseIndices(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  for(int ii=0;ii<ne;ii++)
    {
      Vertex2d v = GetMesh()->vertex2d(indices[ii]);
      T(0,ii) = v.x();               
      T(1,ii) = v.y();
    }
}

/*---------------------------------------------------*/

void DwrFem2d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<2>* I = dynamic_cast<const IntegratorQ1Q2<2>*>(GetIntegrator());
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
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/*---------------------------------------------------*/

void DwrFem2d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<2>* I = dynamic_cast<const IntegratorQ1Q2<2>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation  (TH,iq);
      TransformationQ1(TL,iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,iq);
      I->AdjointForm(EQ,__F,HighOrderFem,LowOrderFem,__U,__Q);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFem2d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  const IntegratorQ1Q2<2>* I = dynamic_cast<const IntegratorQ1Q2<2>*>(GetIntegrator());
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
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}
}
