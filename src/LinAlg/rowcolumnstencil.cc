#include  "rowcolumnstencil.h"
#include  "sparsestructure.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{
void RowColumnStencil::memory(int n, int nt)
{
  srow.reservesize(n);
  ColumnStencil::memory(n,nt);
}

/*-------------------------------------------------------------*/

void RowColumnStencil::memory(const SparseStructureInterface* SI)
{
  ColumnStencil::memory(SI);

  const SparseStructure* SS = dynamic_cast<const SparseStructure*>(SI);
  assert(SS);

  for(int i=0;i<SS->n();i++)
    {
      row(i) = i;
    }
}
}
