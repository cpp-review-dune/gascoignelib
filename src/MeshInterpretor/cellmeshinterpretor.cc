#include  "cellmeshinterpretor.h"
#include  "visudatacompvector.h"
#include  "sparsestructure.h"

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

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const std::vector<Vertex2d>& p0, int comp, const nvector<double>& d) const
{
  int count = 0;
  for (int i=0; i<p0.size(); i++)
    {
      count += RhsPoint(f,p0[i],comp,d[i]);
    }
  return count;
}

/* ----------------------------------------- */

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const std::vector<Vertex3d>& p0, int comp, const nvector<double>& d) const
{
  int count = 0;
  for (int i=0; i<p0.size(); i++)
    {
      count += RhsPoint(f,p0[i],comp,d[i]);
    }
  return count;
}

/*-----------------------------------------*/

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Vertex2d& p0, int comp, double d) const
{
  assert(0);
//   nvector<Vertex2d> p(4);
//   Vertex2d Tranfo_p0;
   
//   int iq=cell_number(p0,p);
//   if (iq==-1) return 0;
  
//   Trafo(p,p0,Tranfo_p0);
  
//   Transformation(T,iq);
//   GetFem()->ReInit(T);
//   GetIntegrator()->RhsPoint(F,*GetFem(),Tranfo_p0,comp);

//   LocalToGlobal(f,iq,__F,d);

//   return 1;
}

/* ----------------------------------------- */

int CellMeshInterpretor::RhsPoint(GlobalVector& f, const Vertex3d& p0, int comp, double d) const
{
  assert(0);
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
