#include  "pointmatrix.h"
#include  "simplesparsestructureadaptor.h"
#include  "nodesparsestructureadaptor.h"
#include  "componentsparsestructureadaptor.h"
#include  "giota.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
PointMatrix::PointMatrix(int ncomp, string type) : MatrixInterface(), _ncomp(ncomp)  
{
  if(type=="node")
    {
//       SSAP = new SimpleSparseStructureAdaptor;
      SSAP = new NodeSparseStructureAdaptor(_ncomp);
    }
  else if(type=="component")
    {
      SSAP = new ComponentSparseStructureAdaptor(_ncomp);
    }
  else
    {
      cerr << "PointMatrix::PointMatrix(): unknown type "<< type<<endl;
      abort();
    }

}

/* ----------------------------------------- */

PointMatrix::~PointMatrix() 
{
  if(SSAP) {delete SSAP; SSAP=NULL;}
}

/* ----------------------------------------- */

void PointMatrix::ReInit(const SparseStructureInterface* S)
{
  SSAP->InitStructure(S);
  SimpleMatrix::ReInit(SSAP->n(),SSAP->nentries());
  SSAP->FillStencil(ST);
}

/* ----------------------------------------- */

void PointMatrix::vmult(GlobalVector& y, const GlobalVector& x, double d) const
{
  assert(SSAP->GetName()=="Node");
  SimpleMatrix::vmult(y,x,d);
}

/* ----------------------------------------- */

void PointMatrix::vmult_transpose(GlobalVector& y, const GlobalVector& x, double d) const
{
  assert(SSAP->GetName()=="Node");
  SimpleMatrix::vmult_transpose(y,x,d);
}

/*-----------------------------------------*/

void PointMatrix::dirichlet(int inode, const vector<int>& cv)
{
  assert(SSAP);
  SimpleMatrix::dirichlet(SSAP->GetIndicesDirichlet(inode,cv));
}

/*-----------------------------------------*/

void PointMatrix::dirichlet_only_row(int inode, const vector<int>& cv)
{
  assert(SSAP);
  SimpleMatrix::dirichlet_only_row(SSAP->GetIndicesDirichlet(inode,cv));
}

/*-----------------------------------------*/

void PointMatrix::entry_diag(int i, const nmatrix<double>& M)
{
  IntVector cv(_ncomp); iota(cv.begin(),cv.end(),0);
  dirichlet(i,cv);
}

/*-----------------------------------------*/

void PointMatrix::entry(const IntVector::const_iterator start, const IntVector::const_iterator stop, const EntryMatrix& M, double s)
{
  int n = stop-start;

  for(int ii=0;ii<n;ii++)
    {
      int i = *(start+ii);
      for(int c=0;c<_ncomp;c++)
	{
	  int iglob = SSAP->index(i,c);
	  for(int jj=0;jj<n;jj++)
	    {
	      int j = *(start+jj);
	      for(int d=0;d<_ncomp;d++)
		{
		  int jglob = SSAP->index(j,d);
		  int pos = ST.Find(iglob,jglob);
		  value[pos] += s*M(ii,jj,c,d);
		}
	    }
	}
    }
}

/*-----------------------------------------*/

void PointMatrix::AddMassWithDifferentStencil(const MatrixInterface* MP, const TimePattern& TP, double s)
{
  const SimpleMatrix* SM = dynamic_cast<const SimpleMatrix*>(MP);
  assert(SM);

  int n = SSAP->nnodes();

  const ColumnStencil*  SMS = dynamic_cast<const ColumnStencil*>(SM->GetStencil());
  const NodeSparseStructureAdaptor* NSMS = dynamic_cast<const NodeSparseStructureAdaptor*>(SSAP);
  assert(NSMS);

  assert(n=SMS->n());

  for(int i=0;i<n;i++)
    {
      for(int pos=SMS->start(i);pos<SMS->stop(i);pos++)
	{
	  int j = SMS->col(pos);
	  double m = s * SM->GetValue(pos);

	  for(int c=0;c<_ncomp;c++)
	    {
	      int isystem = NSMS->index(i, c);
	      for(int d=0;d<_ncomp;d++)
		{
		  int jsystem = NSMS->index(j, d);
		  bool found=0;
		  for(int pos2=ST.start(isystem);pos2<ST.stop(isystem);pos2++)
		    {
		      if(ST.col(pos2)==jsystem)
			{
			  found = 1;
			  value[pos2] += m * TP(c,d);
			}
		    }
		  if(!found) cerr << "not found ";
		}
	    }
	}
    }


//   for(int i=0;i<n;i++)
//     {
//       for(int c=0;c<_ncomp;c++)
// 	{
// 	  intw iglob = SSAP->index(i,c);
// 	  for(int j=0;j<n;j++)
// 	    {
// 	      for(int d=0;d<_ncomp;d++)
// 		{
// 		  int jglob = SSAP->index(j,d);
// 		  int pos = ST.Find(iglob,jglob);
// 		  value[pos] += s * TP(c,d) * SM->GetValue(i,j);
// 		}
// 	    }
// 	}
//     }
}

/*-----------------------------------------*/

void PointMatrix::RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah)
{
  assert(0);

//   const UnstructuredStencil& USH = GetStencil();
//   const UnstructuredStencil& USh = Ah.GetStencil();

//   const UnstructuredStencil& us = I.GetStencil();

//   for(int k=0;k<USh.n();k++)
//     {
//       for(int pos2=us.start(k);pos2<us.stop(k);pos2++)
// 	{	      
// 	  int i = us.col(pos2);
// 	  double alpha_ik = I.Alpha(pos2);

// 	  for(int pos=USh.start(k);pos<USh.stop(k);pos++)
// 	    {
// 	      int l = USh.col(pos);
	      
// 	      for(int pos3=us.start(l);pos3<us.stop(l);pos3++)
// 		{
// 		  int j = us.col(pos3);
// 		  double alpha_jl = I.Alpha(pos3);
// 		  int pos4 = USH.Find(i,j);
// 		  double d = alpha_ik*alpha_jl;

// 		  mat(pos4)->add(d, *Ah.mat(pos));
// 		}
// 	    }
// 	}
//     }
}
}
