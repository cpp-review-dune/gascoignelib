#include  "cellmeshinterpretor.h"
#include  "visudatacompvector.h"
#include  "sparsestructure.h"
#include  "pointfunctional.h"

using namespace std;

/* ----------------------------------------- */

void CellMeshInterpretor::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetMesh()->ncells();iq++)
    {
      nvector<int> indices = GetLocalIndices(iq);
      S->build_add(indices.begin(), indices.end());
    }
  S->build_end();  
}

/* ----------------------------------------- */

void CellMeshInterpretor::Transformation(FemInterface::Matrix& T, int iq) const
{
  int dim = GetMesh()->dimension();
  int ne = GetMesh()->nodes_per_cell();

  nvector<int> indices = GetMesh()->IndicesOfCell(iq);
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

  CompVector<double> lerr(ncomp,3); 

//   int dim = GetMesh()->dimension();
//   int ne = GetMesh()->nodes_per_cell();
//   nmatrix<double> T(dim, ne);

  nmatrix<double> T;
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
	  err(2,c) = GascoigneMath::max(err(2,c),lerr(2,c));
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

void CellMeshInterpretor::RhsNeumann(GlobalVector& f, const Equation& EQ, const IntSet& Colors,  const NeumannData& NRHS, double s) const
{
  nmatrix<double> T;
  for(IntSet::const_iterator p=Colors.begin();p!=Colors.end();p++)
    {
      int col = *p;
      const IntVector& q = GetMesh()->CellOnBoundary(col);
      const IntVector& l = GetMesh()->LocalOnBoundary(col);
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

double CellMeshInterpretor::compute_element_mean_matrix(int iq, EntryMatrix& E) const
{
  assert(0);
}

/* ----------------------------------------- */

double CellMeshInterpretor::PressureFilter(nvector<double>& PF) const
{
  int nv = GetMesh()->nodes_per_cell(); // = 4 oder 8
  EntryMatrix  E(nv,1);

  PF.resize(n());
  PF.zero();
  double omega = 0.;
  for(int iq=0;iq<GetMesh()->ncells();++iq)
    {
      omega += compute_element_mean_matrix(iq,__E);

      nvector<int> ind = GetMesh()->IndicesOfCell(iq);

      for(int i=0;i<ind.size();i++)
	{
	  for(int j=0;j<ind.size();j++)
	    {
	      PF[ind[j]] += __E(i,j,0,0);
	    }
	}
    }
  return omega;
}

/* ----------------------------------------- */

void CellMeshInterpretor::DiracRhs(GlobalVector& f, const RightHandSideData& RHS, double s) const
{
  assert(0);
}

/* ----------------------------------------- */

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Functional* F) const
{
  const PointFunctional* FP = dynamic_cast<const PointFunctional*>(F);
  assert(FP);

  int comp = FP->GetComp();
  const nvector<double>& d = FP->weights();
  int count = 0;

  if (GetMesh()->dimension()==2)
    {
      const std::vector<Vertex2d>& p0 = FP->points2d();
      
      for (int i=0; i<p0.size(); i++)
	{
	  count += RhsPoint(f,p0[i],comp,d[i]);
	}
    }
  else
    {
      const std::vector<Vertex3d>& p0 = FP->points3d();
      
      for (int i=0; i<p0.size(); i++)
	{
	  count += RhsPoint(f,p0[i],comp,d[i]);
	}
    }
  return count;
}

/*-----------------------------------------*/

void CellMeshInterpretor::VertexTransformation(const nvector<Vertex2d>& p, 
					       const Vertex2d& p0, Vertex2d& tp) const
{
  std::vector<double> DT(4);

  double x0 = p[0].x();
  double x1 = p[1].x();
  double x2 = p[2].x();
  double x3 = p[3].x();

  double y0 = p[0].y();
  double y1 = p[1].y();
  double y2 = p[2].y();
  double y3 = p[3].y();

  double a1,a2,b1,b2,c1,c2,d1,d2;
  
  
  d1 = x0; a1 = x1-x0; b1 = x3-x0; c1 = x2-x1+x0-x3;
  d2 = y0; b2 = y3-y0; a2 = y1-y0; c2 = y2-y1-y3+y0;
  
  tp.x()=0.5;
  tp.y()=0.5;
  
  Vertex2d dp;
  Vertex2d res;

  res.x() = p0.x() - ( a1*tp.x() + b1*tp.y() + c1*tp.x()*tp.y() + d1);
  res.y() = p0.y() - ( a2*tp.x() + b2*tp.y() + c2*tp.x()*tp.y() + d2);
  
  do 
    {
      DT[0] = a1+c1*tp.y();
      DT[1] = b1+c1*tp.x();

      DT[2] = a2+c2*tp.y();
      DT[3] = b2+c2*tp.x();

      double det = DT[0]*DT[3]-DT[1]*DT[2];
      assert(fabs(det)>=1.e-20);
      
      dp.x()=(DT[3]*res.x()-DT[2]*res.y())/det;
      dp.y()=(DT[0]*res.y()-DT[1]*res.x())/det;

      tp.x()+=dp.x();
      tp.y()+=dp.y();

      res.x() = p0.x() - ( a1*tp.x() + b1*tp.y() + c1*tp.x()*tp.y() + d1);
      res.y() = p0.y() - ( a2*tp.x() + b2*tp.y() + c2*tp.x()*tp.y() + d2);
      assert(tp.norm_l8()<=2);
    } while (res.norm_l8()>1.e-14);
}

/*-----------------------------------------*/

int CellMeshInterpretor::GetCellNumber(const Vertex2d& p0, nvector<Vertex2d>& p) const
{
  nvector<double> a0(2),a1(2),b(2),n(2);
  
  int nv = GetMesh()->nodes_per_cell(); // = 4 oder 8
  p.resize(nv);
  for(int iq=0; iq<GetMesh()->ncells(); ++iq)
    {
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

  nvector<Vertex2d> p(GetMesh()->nodes_per_cell());
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
}

/* ----------------------------------------- */

double CellMeshInterpretor::ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const 
{
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
