#include  "simplematrix.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
ostream& SimpleMatrix::Write(ostream& os) const
{
  int n = ST.n();
  for(int i=0;i<n;i++)
    {
      os << i << endl;
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
        {
          os << ST.col(pos) << ": " << value[pos] << ", ";
        }
      os << endl;
    }

  return os;
}

/* ----------------------------------------- */

void SimpleMatrix::ReInit(int n, int nentries)
{
  ST.memory(n,nentries);
  value.reservesize(nentries);
}

/* ----------------------------------------- */

void SimpleMatrix::ReInit(const SparseStructureInterface* SI)
{
  const SparseStructure* S = dynamic_cast<const SparseStructure*>(SI);
  SimpleMatrix::ReInit(S->n(),S->ntotal());
  ST.memory(S);
}

/* ----------------------------------------- */

void SimpleMatrix::entry(niiterator start, niiterator stop, const EntryMatrix& M, double s)
{
  int n = stop-start;

  for(int ii=0;ii<n;ii++)
    {
      int i = *(start+ii);
      for(int jj=0;jj<n;jj++)
        {
          int j = *(start+jj);
          int pos = ST.Find(i,j);
          value[pos] += s*M(ii,jj,0,0);
        }
    }
}

/* ----------------------------------------- */

void SimpleMatrix::vmult_time(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s) const
{
  int n = ST.n();
  assert(n==y.n());
  assert(n==x.n());
  assert(x.ncomp()==y.ncomp());

  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
        {
          int j = ST.col(pos);
          for(int c=0;c<x.ncomp();c++)
            {
              for(int d=0;d<x.ncomp();d++)
                {
                  y(i,c) += s*value[pos]* TP(c,d) * x(j,d);
                }
            }
        }
    }
}

/* ----------------------------------------- */

void SimpleMatrix::vmult(DoubleVector& y, const DoubleVector& x, double d) const
{
  int n = ST.n();
  assert(n==y.size());
  assert(n==x.size());

  DoubleVector::iterator py=y.begin();
  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
        {
          int j = ST.col(pos);
          // y[i] += d*value[pos]*x[j];
          *py += d*value[pos]*x[j];
        }
      py++;
    }
}

/* ----------------------------------------- */

void SimpleMatrix::vmult_transpose(DoubleVector& y, const DoubleVector& x, double d) const
{
  int n = ST.n();
  assert(n==y.size());
  assert(n==x.size());

  DoubleVector::const_iterator px=x.begin();
  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
// 	  y[j] += d*value[pos]*x[i];
	  y[j] += d*value[pos] * *px;
	}
      px++;
    }
}

/*-----------------------------------------*/

void SimpleMatrix::vmult_comp(int c, int d, GlobalVector& y, const GlobalVector& x, double s) const
{
  int n = ST.n();
  assert(n==y.n());
  assert(n==x.n());

  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  y(i,c) += s*value[pos]*x(j,d);
	}
    }
}

/*-----------------------------------------*/

void SimpleMatrix::vmult_comp_trans(int c, int d, GlobalVector& y, const GlobalVector& x, double s) const
{
  int n = ST.n();
  assert(n==y.n());
  assert(n==x.n());

  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  y(j,c) += s*value[pos]*x(i,d);
	}
    }
}

/*-----------------------------------------*/

void SimpleMatrix::dirichlet(const IntVector& indices)
{
  for(int ii=0;ii<indices.size();ii++)
    {
      int i = indices[ii];
      if(i<0) cerr << "SimpleMatrix::dirichlet indices: " << indices << endl;
      assert(i>=0);

      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  if(j==i) value[pos]=1.;
	  else 
	    {
	      value[pos] = 0.;
	      for(int pos2=ST.start(j);pos2<ST.stop(j);pos2++)
		{
		  if(ST.col(pos2)==i) 
		    {
		      value[pos2]=0.;
		      break;
		    }
		}
	    }
	}
    }
}

/*-----------------------------------------*/

void SimpleMatrix::dirichlet_only_row(const IntVector& indices)
{
  for(int ii=0;ii<indices.size();ii++)
    {
      int i = indices[ii];
      if(i<0) cerr << "SimpleMatrix::dirichlet indices: " << indices << endl;
      assert(i>=0);

      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  if(j==i) value[pos]=1.;
	  else 
	    {
	      value[pos] = 0.;
	    }
	}
    }
}

/*-----------------------------------------*/

void SimpleMatrix::transpose()
{
  for(int i=0; i<ST.n(); i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  if(j<i)
	    {
	      double help = value[pos];
	      for(int pos2=ST.start(j);pos2<ST.stop(j);pos2++)
		{
		  if(ST.col(pos2)==i)
		  {
		    value[pos] = value[pos2];
		    value[pos2] = help;
		  }
		}
	    }
	}
    }
}

/*-----------------------------------------*/

void SimpleMatrix::entry_diag(int i, const nmatrix<double>& M)
{
  int pos = ST.Find(i,i);
  value[pos] = M(0,0);
}

/*-----------------------------------------*/

void SimpleMatrix::PrepareJacobi(double s)
{
  int n = ST.n();
  
  _diag.resize(n);
  
  for(int i=0; i<n; i++)
  {
    _diag[i] = s * value[ST.Find(i,i)];
  }
}

/*-----------------------------------------*/

void SimpleMatrix::JacobiVector(GlobalVector &y) const
{
  int n = ST.n();
  assert(n==y.n());
  
  for(int i=0; i<n; i++)
  {
    for(int c=0; c<y.ncomp(); c++)
    {
      y(i,c) *= sqrt(_diag[i]);
    }
  }
}

/*-----------------------------------------*/

void SimpleMatrix::Jacobi(GlobalVector &y) const
{
  int n = ST.n();
  assert(n==y.n());

  for(int i=0; i<n; i++)
  {
    for(int c=0; c<y.ncomp(); c++)
    {
      y(i,c) /= sqrt(_diag[i]);
    }
  }
}

/*-----------------------------------------*/

void SimpleMatrix::vmult_time_Jacobi(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s) const
{
  int n = ST.n();
  assert(n==y.n());
  assert(n==x.n());
  assert(x.ncomp()==y.ncomp());

  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
	  for(int c=0;c<x.ncomp();c++)
	    {
	      for(int d=0;d<x.ncomp();d++)
		{
		  y(i,c) += s * value[pos]* TP(c,d) * x(j,d) / sqrt(_diag[i] * _diag[j]);
		}
	    }
	}
    }
}

}
