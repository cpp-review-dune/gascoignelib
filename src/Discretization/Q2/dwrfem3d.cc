#include "dwrfem.h" 
#include "galerkinintegratorq2.h"
#include "baseq23d.h"
#include "integratorq1q2.h"
#include <ext/hash_set>

using namespace std;

namespace Gascoigne
{

/*---------------------------------------------------*/

DwrFem3d::DwrFem3d() : Q23d()
{
}

/*---------------------------------------------------*/

void DwrFem3d::BasicInit(const ParamFile*  paramfile)
{
  assert(PatchDiscretization::GetIntegrator()==NULL);
  PatchDiscretization::GetIntegratorPointer() =  new IntegratorQ1Q2<3>;

  GetIntegratorPointer()->BasicInit();

  assert(PatchDiscretization::GetFem()==NULL);
  typedef Transformation3d<BaseQ23d>           TransQ2;
  typedef FiniteElement<3,2,TransQ2,BaseQ23d>  FiniteElement;

  PatchDiscretization::GetFemPointer() =  new FiniteElement;
  PatchDiscretization::BasicInit(paramfile);
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

void DwrFem3d::DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  Vertex3d Tranfo_p0;
   
  int iq = GetPatchNumber(p0,Tranfo_p0);
  if (iq==-1)
    {
      cerr << "DwrFem3d::DiracRhsPoint point not found\n";
      abort();
    }

  nmatrix<double> TH,TL;

  Transformation  (TH,iq);
  TransformationQ1(TL,iq);

  const FemInterface& HighOrderFem(*GetFem());

  HighOrderFem.ReInit(TH);
  LowOrderFem .ReInit(TL);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);
  
  GlobalToLocalData(iq);
  GlobalToGlobalData();
  DRHS.SetParameterData(__qq);

  I->DiracRhsPoint(__F,HighOrderFem,LowOrderFem,Tranfo_p0,DRHS,i,__Q);
  PatchDiscretization::LocalToGlobal(f,__F,iq,s);
}

/*---------------------------------------------------*/

void DwrFem3d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__qq);

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
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFem3d::AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  EQ.SetParameterData(__qq);

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
      I->AdjointForm(EQ,__F,HighOrderFem,LowOrderFem,__U,__Q);
      PatchDiscretization::LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void DwrFem3d::BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, 
    const BoundaryEquation& BE, double d) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  BE.SetParameterData(__qq);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<nvector<int> >& patch2cell  =
      GetGascoigneMesh()->GetPatchIndexHandler().GetAllPatch2Cell();

  nvector<int> cell2patch(GetMesh()->ncells());
  for (int p=0;p<patch2cell.size();++p)
    for (int i=0;i<patch2cell[p].size();++i)
       cell2patch[patch2cell[p][i]]=p;
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;

    __gnu_cxx::hash_set<int> habschon;
    
    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2patch[iq];
      
      // gabs den patch schon?
      if (habschon.find(ip)!=habschon.end()) continue;
      habschon.insert(ip);

      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocal(__U,u,ip);
      I->BoundaryForm(BE,__F,HighOrderFem,LowOrderFem,__U,ile,col,__Q);
      PatchDiscretization::LocalToGlobal(f,__F,ip,d);
    }
  }
}

/* ----------------------------------------- */

void DwrFem3d::Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  RHS.SetParameterData(__qq);

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
      PatchDiscretization::LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void DwrFem3d::BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
{
  nmatrix<double> TH,TL;

  GlobalToGlobalData();
  NRHS.SetParameterData(__qq);
  
  const IntegratorQ1Q2<3>* I = dynamic_cast<const IntegratorQ1Q2<3>*>(GetIntegrator());
  assert(I);

  const FemInterface& HighOrderFem(*GetFem());

  const nvector<nvector<int> >& patch2cell  =
      GetGascoigneMesh()->GetPatchIndexHandler().GetAllPatch2Cell();

  nvector<int> cell2patch(GetMesh()->ncells());
  for (int p=0;p<patch2cell.size();++p)
    for (int i=0;i<patch2cell[p].size();++i)
       cell2patch[patch2cell[p][i]]=p;
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
  {
    int col = *p;
    __gnu_cxx::hash_set<int> habschon;

    const IntVector& q = *GetMesh()->CellOnBoundary(col);
    const IntVector& l = *GetMesh()->LocalOnBoundary(col);
    for (int i=0; i<q.size(); i++)
    {
      int iq  = q[i];
      int ip  = cell2patch[iq];
      
      // gabs den patch schon?
      if (habschon.find(ip)!=habschon.end()) continue;
      habschon.insert(ip);
      
      int ile = l[i];

      Transformation  (TH,ip);
      TransformationQ1(TL,ip);

      HighOrderFem.ReInit(TH);
      LowOrderFem .ReInit(TL);

      GlobalToLocalData(ip);
      I->BoundaryRhs(NRHS,__F,HighOrderFem,LowOrderFem,ile,col,__Q);
      PatchDiscretization::LocalToGlobal(f,__F,ip,s);
    }
  }
}

/* ----------------------------------------- */

void Gascoigne::DwrFem3d::MassMatrix(MatrixInterface& A) const
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

      I->MassMatrix(__E,HighOrderFem,LowOrderFem);
      PatchDiscretization::LocalToGlobal(A,__E,iq,1.);
    }
}
}
