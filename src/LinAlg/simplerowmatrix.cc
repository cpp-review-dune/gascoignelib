#include  "simplerowmatrix.h"
#include  "sparsestructure.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
std::ostream& SimpleRowMatrix::Write(std::ostream& os) const
{
  int n = ST.n();
  for(int i=0;i<n;i++)
    {
      os << i << endl;
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  os << value[pos] << " ";
	}
      os << endl;
    }
  return os;
}

/* ----------------------------------------- */

void SimpleRowMatrix::vmult(nvector<double>& y, const nvector<double>& x, double d) const
{
  int n = ST.n();
  assert(n==y.size());
  assert(n==x.size());

  nvector<double>::iterator py=y.begin();
  for(int i=0;i<n;i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  int j = ST.col(pos);
// 	  y[i] += d*value[pos]*x[j];
	  *py += d*value[pos]*x[j];
	}
      py++;
    }
}

/* ----------------------------------------- */

void SimpleRowMatrix::ReInit(int n, int nentries)
{
  ST.memory(n,nentries);
  value.reservesize(nentries);
}

/* ----------------------------------------- */

void SimpleRowMatrix::ReInit(const SparseStructureInterface* SI)
{
  const SparseStructure* S = dynamic_cast<const SparseStructure*>(SI);
  SimpleRowMatrix::ReInit(S->n(),S->ntotal());
  ST.memory(S);
}

/* ----------------------------------------- */

void SimpleRowMatrix::entry(niiterator start, niiterator stop, const EntryMatrix& M, double s)
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
}
