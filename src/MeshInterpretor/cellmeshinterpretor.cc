#include  "cellmeshinterpretor.h"
#include  "visudatacompvector.h"
#include  "sparsestructure.h"
#include  "pointfunctional.h"
#include  "pressurefilter.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
void CellMeshInterpretor::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetMesh()->ncells();iq++)
    {
      IntVector indices = GetLocalIndices(iq);
      S->build_add(indices.begin(), indices.end());
    }
  S->build_end();  
}

/* ----------------------------------------- */

void CellMeshInterpretor::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  IntVector indices = GetMesh()->IndicesOfCell(iq);
  assert(ne==indices.size());

  T.memory(dim,ne);
  if(dim==2)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	}
    }
  else if(dim==3)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex3d v = GetMesh()->vertex3d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	  T(2,ii) = v.z();
	}
    }
}

/* ----------------------------------------- */
/* ----------------------------------------- */
/* ----------------------------------------- */

void CellMeshInterpretor::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  EQ.SetParameterData(__qq);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      
      GetIntegrator()->Form(EQ,__F,*GetFem(),__U,__Q);
      LocalToGlobal(f,__F,iq,d);
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  EQ.SetParameterData(__qq);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      GetIntegrator()->Matrix(EQ,__E,*GetFem(),__U,__Q);
      LocalToGlobal(A,__E,iq,d);
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::MassMatrix(MatrixInterface& A) const
{
  nmatrix<double> T;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GetIntegrator()->MassMatrix(__E,*GetFem());
      LocalToGlobal(A,__E,iq,1.);
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const
{
//   const IntegrationFormulaInterface& IF = ErrorFormula();

  int ncomp = u.ncomp();
  err.ncomp() = ncomp;
  err.reservesize(3);
  err = 0.;

  GlobalVector lerr(ncomp,3); 

  nmatrix<double> T;
  
  GlobalToGlobalData();
  ES->SetParameterData(__qq);
  
  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);
      GlobalToLocal(__U,u,iq);
      GetIntegrator()->ErrorsByExactSolution(lerr,*GetFem(),*ES,__U,__Q);
      for(int c=0;c<ncomp;c++)  
	{
	  err(0,c) += lerr(0,c);
	  err(1,c) += lerr(1,c);
	  err(2,c) = Gascoigne::max(err(2,c),lerr(2,c));
	}
    }
  for(int c=0;c<ncomp;c++)  
    {
      err(0,c) = sqrt(err(0,c));
      err(1,c) = sqrt(err(1,c));
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::Rhs(GlobalVector& f, const RightHandSideData& RHS, double s) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  RHS.SetParameterData(__qq);
  
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocalData(iq);
      GetIntegrator()->Rhs(RHS,__F,*GetFem(),__Q);
      LocalToGlobal(f,__F,iq,s);
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::RhsNeumann(GlobalVector& f, const IntSet& Colors,  const NeumannData& NRHS, double s) const
{
  nmatrix<double> T;
  
  GlobalToGlobalData();
  NRHS.SetParameterData(__qq);
  
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = *GetMesh()->CellOnBoundary(col);
      const IntVector& l = *GetMesh()->LocalOnBoundary(col);
      for (int i=0; i<q.size(); i++)
	{
	  int iq  = q[i];
	  int ile = l[i];

	  Transformation(T,iq);
	  GetFem()->ReInit(T);

	  GlobalToLocalData(iq);
	  GetIntegrator()->RhsNeumann(NRHS,__F,*GetFem(),ile,col,__Q);
	  LocalToGlobal(f,__F,iq,s);
	}
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::InitFilter(DoubleVector& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;
  for(int iq=0; iq<GetMesh()->ncells(); ++iq)
    {
      int nv = GetMesh()->nodes_per_cell(iq);
      EntryMatrix  E(nv,1);

      Transformation(T,iq);
      GetFem()->ReInit(T);

      double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
      PF->AddDomainPiece(cellsize);

      IntVector ind = GetMesh()->IndicesOfCell(iq);

      for(int i=0;i<ind.size();i++)
 	{
	  for(int j=0;j<ind.size();j++)
	    {
	      F[ind[j]] += E(i,j,0,0);
	    }
      	}
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const
{
  const NewDiracRightHandSide* DRHS = dynamic_cast<const NewDiracRightHandSide*>(&RHS);
  assert(DRHS);

  int dim = GetMesh()->dimension();
  vector<int> comps = DRHS->GetComps();
  int nn = comps.size();

  vector<double> up(nn,0);
 
  GlobalToGlobalData();

  if (dim == 2)
    {
      vector<Vertex2d> v2d = DRHS->GetPoints2d();
      assert(nn==v2d.size());
      
      for(int i=0;i<nn;++i)
	{
	  DiracRhsPoint(f,DRHS,v2d[i],i,s);
	}
    }
  else if (dim == 3)
    {
      vector<Vertex3d> v3d = DRHS->GetPoints3d();
      assert(nn==v3d.size());
      for(int i=0;i<nn;++i)
	{
	  DiracRhsPoint(f,DRHS,v3d[i],i,s);
	}
    }
  else
    {
      cerr << "wrong dim = " << dim << endl;
      abort();
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::DiracRhsPoint(GlobalVector& f,const NewDiracRightHandSide* DRHS,const Vertex2d& p0,int i,double s) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  nvector<Vertex2d> p(8);
  Vertex2d Tranfo_p0;
   
  int iq = GetCellNumber(p0,p);
  if (iq==-1)
    {
      cerr << "CellMeshInterpretor::DiracRhsPoint point not found\n";
      abort();
    }

  VertexTransformation(p,p0,Tranfo_p0);

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GetIntegrator()->DiracRhsPoint(__F,*GetFem(),Tranfo_p0,DRHS,i,__Q);
  LocalToGlobal(f,__F,iq,s);
}

/* ----------------------------------------- */

void CellMeshInterpretor::DiracRhsPoint(GlobalVector& f,const NewDiracRightHandSide* DRHS,const Vertex3d& p0,int i,double s) const
{
  cerr << "CellMeshInterpretor::DiracRhsPoint not written in 3d\n";
  abort();
}

/* ----------------------------------------- */

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Functional* F) const
{
  const PointFunctional* FP = dynamic_cast<const PointFunctional*>(F);
  assert(FP);

  int comp = FP->GetComp();
  const DoubleVector& d = FP->weights();
  int count = 0;

  if (GetMesh()->dimension()==2)
    {
      const vector<Vertex2d>& p0 = FP->points2d();
      
      for (int i=0; i<p0.size(); i++)
	{
	  count += RhsPoint(f,p0[i],comp,d[i]);
	}
    }
  else
    {
      const vector<Vertex3d>& p0 = FP->points3d();
      
      for (int i=0; i<p0.size(); i++)
	{
	  count += RhsPoint(f,p0[i],comp,d[i]);
	}
    }
  return count;
}

/*-----------------------------------------*/

int CellMeshInterpretor::GetCellNumber(const Vertex2d& p0, nvector<Vertex2d>& p) const
{
  DoubleVector a0(2),a1(2),b(2),n(2);
  
  for(int iq=0; iq<GetMesh()->ncells(); ++iq)
    {
      int nv = GetMesh()->nodes_per_cell(iq); // = 4 oder 8
      p.resize(nv);
      for (int j=0; j<nv; ++j)
	{    
	  p[j]=GetMesh()->vertex2d(GetMesh()->vertex_of_cell(iq,j));
	}
      int s=0;
      for (int j=0;j<nv;++j)    
	{
	  n[0] =-p[(j+1)%4].y()+p[j].y();
	  n[1] = p[(j+1)%4].x()-p[j].x();
	  b[0] = p0.x()-p[j].x();
	  b[1] = p0.y()-p[j].y();
	  if (n*b>=-1.e-15) ++s;
	}
      if (s==nv) return iq;
    }
  return -1;
}

/*-----------------------------------------*/

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Vertex2d& p0, int comp, double d) const
{
  __F.ReInit(f.ncomp(),GetFem()->n());

  nvector<Vertex2d> p(8);
  Vertex2d Tranfo_p0;
   

  int iq = GetCellNumber(p0,p);
  if (iq==-1) return 0;

  VertexTransformation(p,p0,Tranfo_p0);

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GetIntegrator()->RhsPoint(__F,*GetFem(),Tranfo_p0,comp);
  LocalToGlobal(f,__F,iq,d);

  return 1;
}

/* ----------------------------------------- */

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Vertex3d& p0, int comp, double d) const
{
  // p0 must be mesh point !!!!
  int found = -1;
  for (int i=0; i<GetMesh()->nnodes(); i++)
    {
      if (GetMesh()->vertex3d(i)==p0) found = i;
    }
  if (found<0)
  {
    cerr << "point is not a mesh point!!!\n";
    return 0;
  }

  f(found,comp) = d;
  
  return 1;
}

/* ----------------------------------------- */

double CellMeshInterpretor::ComputeBoundaryFunctional(const GlobalVector& u, const BoundaryFunctional& BF) const 
{
  assert(0);
  return 0.;
}

/* ----------------------------------------- */

double CellMeshInterpretor::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const 
{
  GlobalToGlobalData();
  F.SetParameterData(__qq);
  
  nmatrix<double> T;
  double j=0.;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      GlobalToLocal(__U,u,iq);
      j += GetIntegrator()->ComputeDomainFunctional(F,*GetFem(),__U,__Q);
    }
  return j;
}

/* ----------------------------------------- */

double CellMeshInterpretor::ComputeNewPointFunctional(const GlobalVector& u, const NewPointFunctional* FP) const
{
  int dim = GetMesh()->dimension();
  vector<int> comps = FP->GetComps();
  int nn = comps.size();

  vector<double> up(nn,0);
 
  if (dim == 2)
    {
      vector<Vertex2d> v2d = FP->GetPoints2d();
      assert(nn==v2d.size());
      
      for(int i=0;i<nn;++i)
	{
	  up[i] = ComputePointValue(u,v2d[i],comps[i]);
	}
    }
  else if (dim == 3)
    {
      vector<Vertex3d> v3d = FP->GetPoints3d();
      assert(nn==v3d.size());
      for(int i=0;i<nn;++i)
	{
	  up[i] = ComputePointValue(u,v3d[i],comps[i]);
	}
    }
  else
    abort();

  return FP->J(up);
}
/* ----------------------------------------- */

double CellMeshInterpretor::ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const
{
  nvector<Vertex2d> p(8);
  Vertex2d Tranfo_p0;

  int iq = GetCellNumber(p0,p);
  if (iq==-1)
    {
      cerr << "CellMeshInterpretor::ComputePointValue point not found\n";
      abort();
    }

  VertexTransformation(p,p0,Tranfo_p0);

  nmatrix<double> T;
  Transformation(T,iq);
  GetFem()->ReInit(T);

  GlobalToLocal(__U,u,iq);

  return GetIntegrator()->ComputePointValue(*GetFem(),Tranfo_p0,__U,comp);
}

/* ----------------------------------------- */

double CellMeshInterpretor::ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const
{
  cerr << "ComputePointValue for 3d not written\n";
  assert(0);
}

/* ----------------------------------------- */

/* ----------------------------------------- */

void CellMeshInterpretor::Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell(iq);

  IntVector indices = HM->GetVertices(iq);
  assert(ne==indices.size());
  swapIndices(indices);

  T.memory(dim,ne);
  if(dim==2)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex2d v = GetMesh()->vertex2d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	}
    }
  else if(dim==3)
    {
      for(int ii=0;ii<ne;ii++)
	{
	  Vertex3d v = GetMesh()->vertex3d(indices[ii]);
	  T(0,ii) = v.x();               
	  T(1,ii) = v.y();
	  T(2,ii) = v.z();
	}
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const
{
  IntVector indices = HM->GetVertices(iq);
  swapIndices(indices);

  U.ReInit(u.ncomp(),indices.size());
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      U.equ_node(ii,i,u);
    }
}

/* ----------------------------------------- */

void CellMeshInterpretor::swapIndices(IntVector& indices) const
{
  assert(indices.size()>=4);

  int help = indices[2];
  indices[2] = indices[3];
  indices[3] = help;
  if (indices.size()==8)
    {
      help = indices[6];
      indices[6] = indices[7];
      indices[7] = help;
    }
}
}
