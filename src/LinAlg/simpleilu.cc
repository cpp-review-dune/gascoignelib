#include  "simpleilu.h"


using namespace std;

/* ----------------------------------------- */

void SimpleIlu::ReInit(int n, int nentries)
{
  SimpleMatrix::ReInit(n,nentries);
  p.reservesize(n,-1);
  q.reservesize(n,-1);
}

/*-----------------------------------------*/

void SimpleIlu::solve(nvector<double>& x) const
{
  hin(x);
  forward ();
  backward();
  her(x);
}

/*-----------------------------------------*/

void SimpleIlu::solve_transpose(nvector<double>& x) const
{
  hin(x);
  forward_transpose ();
  backward_transpose();
  her(x);
}

/*-------------------------------------------------------------*/

void SimpleIlu::hin(const nvector<double>& x) const
{
  int n = x.size();
  yp.reservesize(n);
  assert(n==ST.n());
  for(int i=0;i<n;i++)  yp[i] = x[p[i]];
}

/*-------------------------------------------------------------*/

void SimpleIlu::her(nvector<double>& x) const
{
  for(int i=0;i<ST.n();i++)  x[i] = yp[q[i]];
}

/*-------------------------------------------------------------*/

void SimpleIlu::forward() const
{
  for(int i=1; i<ST.n(); i++)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.start(i); pos<ende; pos++)
	{
	  int j = ST.col(pos);
	  yp[i] -= value[pos]*yp[j];
	}
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::backward() const
{
  for(int i=ST.n()-1; i>=0; i--)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.stop(i)-1; pos>ende; pos--)
	{
	  int j = ST.col(pos);
	  yp[i] -= value[pos]*yp[j];
	}
      yp[i]  *= value[ende];
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::forward_transpose() const
{
  for(int i=0; i<ST.n(); i++)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.start(i); pos<ende; pos++)
	{
	  int j = ST.col(pos);
	  int pos2 = ST.Find(j,i);
	  yp[i] -= value[pos2]*yp[j];
	}
      yp[i]  *= value[ende];
    }
}

/*-------------------------------------------------------------*/

void SimpleIlu::backward_transpose() const
{
  for(int i=ST.n()-1; i>=0; i--)
    {
      int ende = ST.diag(i); 
      for(int pos=ST.stop(i)-1; pos>ende; pos--)
	{
	  int j = ST.col(pos);
	  int pos2 = ST.Find(j,i);
	  yp[i] -= value[pos2]*yp[j];
	}
    }
}

/* ----------------------------------------- */
  
void SimpleIlu::compute_ilu()
{
  for(int i=0; i<ST.n(); i++)
    {
      for (int pk=ST.start(i); pk<ST.diag(i); pk++)
	{
	  int k = ST.col(pk);

	  value[pk] *= value[ST.diag(k)];

	  for (int pj=ST.diag(k)+1; pj<ST.stop(k); pj++)
	    {
	      int j  = ST.col(pj);
	      // suche ph
	      for (int ph=ST.start(i); ph<ST.stop(i); ph++)
		{
		  if (ST.col(ph)==j)
		    {
		      value[ph] -= value[pk]*value[pj];
		      break;
		    }
		}
	    }
	}
      double d = value[ST.diag(i)];
      value[ST.diag(i)] = 1./d;
    }
}

/*-------------------------------------------------*/

void SimpleIlu::copy_entries(const MatrixInterface*  A)
{
  const SimpleMatrix* AP = dynamic_cast<const SimpleMatrix*>(A);
  assert(AP);

  const ColumnDiagStencil* AS = dynamic_cast<const ColumnDiagStencil*>(AP->GetStencil());
  assert(AS);

  for(int i=0;i<ST.n();i++)
    {
      int pi = p[i];

      for(int posA=AS->start(pi); posA<AS->stop(pi); posA++)
	{
	  int j   = AS->col(posA);
	  int pj  = q[j];
	  bool found=0;
	  for(int pos=ST.start(i); pos<ST.stop(i); pos++)
	    {
	      int k = ST.col(pos);
	      if(k==pj)	
		{
		  value[pos] += AP->GetValue(posA);
		  found=1;
		  break;
		}
	    }
	  if(!found)
	    {
	      std::cout << "not found " << std::endl;
	      abort();
	    }
	}
    }
}
