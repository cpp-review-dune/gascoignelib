#include  "q1.h"
#include  "gascoignemesh.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

Q1::Q1() : CellMeshInterpretor(), HN(NULL) 
{
}

/* ----------------------------------------- */

Q1::~Q1()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

int Q1::n() const
{
  GetMesh()->nnodes();
}

/* ----------------------------------------- */

void Q1::ReInit(const MeshInterface* MP)
{
  CellMeshInterpretor::ReInit(MP);
  HN->ReInit(MP);
}

/* ----------------------------------------- */

void Q1::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  HN->CondenseHanging(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}

/* ----------------------------------------- */

nvector<int> Q1::GetLocalIndices(int iq) const 
{
  nvector<int> indices = GetMesh()->IndicesOfCell(iq);
  return indices;
}

/* ----------------------------------------- */

void Q1::StrongDirichletMatrix(MatrixInterface& A, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  const IntVector& bv = GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      A.dirichlet(bv[i], comp);
    }  
}

/* ----------------------------------------- */

void Q1::StrongDirichletVectorZero(GlobalVector& u, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  const IntVector& bv = GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];
      for(int iii=0;iii<comp.size();iii++)
	{
	  u( index,comp[iii] ) = 0.;
	}
    }  
}

/* ----------------------------------------- */

void Q1::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
{
  const IntVector& vo2n = GetMesh()->Vertexo2n();
  assert(vo2n.size()==uold.n());

  DoubleVector habschon(GetMesh()->nnodes(),0.);  
  nvector<bool> oldnode(GetMesh()->nnodes(),0);

  u.zero();
  for(int i=0;i<vo2n.size();i++)
    {
      int in = vo2n[i];

      if(in>=0) 
	{
	  u.equ_node(in,1.,i,uold);
	  oldnode[in] = 1;
	}
    }

  nmatrix<double> w = GetLocalInterpolationWeights();

  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      IntVector v = GetMesh()->IndicesOfCell(iq);
      for(int iol=0; iol<v.size(); iol++)
	{
	  int io = v[iol];

	  if (oldnode[io])
	    {
	      for (int inl=0; inl<v.size(); inl++)
		{
		  if (iol==inl)        continue;
		  int in = v[inl];
		  if (oldnode[in])     continue;
		  if (habschon[in]>=1) continue;

		  double weight = w(iol,inl);

		  u.add_node(in,weight,io,uold);

		  habschon[in] += weight;
		}
	      
	    }
	}
    }
}

/* ----------------------------------------- */

void Q1::HNAverage(GlobalVector& x) const
{
  HN->Average(x);
}

/* ----------------------------------------- */

void Q1::HNDistribute(GlobalVector& x) const
{
  HN->Distribute(x);
}

/* ----------------------------------------- */

void Q1::HNZero(GlobalVector& x) const
{
  HN->Zero(x);
}

/* ----------------------------------------- */

bool Q1::HNZeroCheck(const GlobalVector& x) const
{
  return HN->ZeroCheck(x);
}

/* ----------------------------------------- */

void Q1::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetMesh()->ncells();iq++)
    {
      nvector<int> indices = GetLocalIndices(iq);
      HN->CondenseHanging(indices);
      S->build_add(indices.begin(), indices.end());
    }
  HN->SparseStructureDiag(S);
  S->build_end();
}

/* ----------------------------------------- */

void Q1::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  CellMeshInterpretor::Matrix(A,u,EQ,d);

  HN->MatrixDiag(u.ncomp(),A);
}

/* ----------------------------------------- */

void Q1::MassMatrix(MatrixInterface& A) const
{
  CellMeshInterpretor::MassMatrix(A);

  HN->MatrixDiag(1,A);  
}

/* ----------------------------------------- */

double Q1::PressureFilter(nvector<double>& PF) const
{
  assert(0);
}

/* ----------------------------------------- */

