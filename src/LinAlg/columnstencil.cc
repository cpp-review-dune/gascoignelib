#include  "columnstencil.h"
#include  "sparsestructure.h"


using namespace std;

/*-------------------------------------------------------------*/
  
std::ostream& operator<<(std::ostream &s, const ColumnStencil& A)
{
  s << "start:\n"<< A.start() << std::endl;
  s << "col:\n"<< A.col() << std::endl;
  return s;
}

/*-------------------------------------------------------------*/

void ColumnStencil::memory(int n, int nt)
{
  scol  .reservesize(nt);
  sstart.reservesize(n+1);
}

/*-------------------------------------------------------------*/

void ColumnStencil::memory(const SparseStructureInterface* SI)
{
  const SparseStructure* SS = dynamic_cast<const SparseStructure*>(SI);
  assert(SS);

  memory(SS->n(),SS->ntotal());

  sstart[0] = 0;
  for(int i=0;i<SS->n();i++)
    {
      int first = sstart[i];
      sstart[i+1] = first + SS->rowsize(i);
      int id = 0;
      for(std::set<int>::const_iterator p=SS->rowbegin(i);
	  p!=SS->rowend(i);p++)
	{
	  scol[first+id] = *p;
	  id++;
	}
    }
}
