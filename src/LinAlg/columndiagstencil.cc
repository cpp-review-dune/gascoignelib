#include  "columndiagstencil.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{
void ColumnDiagStencil::memory(const SparseStructureInterface* SI)
{
  ColumnStencil::memory(SI);

  sdiag.reservesize(n());

  for(int i=0;i<n();i++)
    {
      for (int pos=sstart[i]; pos<sstart[i+1]; pos++)
        {
          if (scol[pos]==i)
            {
              sdiag[i] = pos;
              break;
            }
        }
    }
}

/*-------------------------------------------------------------*/

void ColumnDiagStencil::memory(int n, int nt)
{
  ColumnStencil::memory(n,nt);
  sdiag.reservesize(n);
}
}
