#include  "twinstencil.h"

using namespace std;

/*-------------------------------------------------------------*/
  
namespace Gascoigne
{
void TwinStencil::diagonalfirst()
{
  // has to be called if later on an ilu factorization is performed

  for(int i=0; i<n(); i++)
    {
      int first = sstart[i];
      int found = -1;
      for (int pos=first; pos<sstart[i+1]; pos++)
	{
	  if (scol[pos]==i)
	    {
	      found=pos;
	      break;
	    }
	}
      if (found==-1)
	{
	  cout << "UnstructuredStencil::diagonal not found " << i << endl;
	  abort();
	}
      for (int pos=found; pos>first; pos--)
	{
	  swap(scol[pos],scol[pos-1]);
	}
    }
}

/*-------------------------------------------------------------*/
  
int TwinStencil::half(int i) const
{
  for (int pos=sstart[i]; pos<sstart[i+1]; pos++)
    {
      if (scol[pos]>i) return pos;
    }
  return sstart[i+1];
}

/*-------------------------------------------------------------*/

void TwinStencil::memory(int n, int nt)
{
  ColumnStencil::memory(n,nt);
}

/*-------------------------------------------------------------*/

void TwinStencil::memory(const SparseStructureInterface* SA)
{
  ColumnStencil::memory(SA);
  diagonalfirst();
}
}
