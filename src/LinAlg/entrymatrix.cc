#include  "entrymatrix.h"

using namespace std;


/*-------------------------------------------------------------*/

namespace Gascoigne
{
ostream& operator<<(ostream& s, const EntryMatrix& A)
{
  s << A.Ndof()  << "\t" << A.Mdof() << endl;
  s << A.Ncomp()  << "\t" << A.Mcomp() << endl;
  for(int i=0;i<A.Ndof();i++)
    {
      for(int j=0;j<A.Mdof();j++)
        {
          s << "\n[" << i << "," << j << "] ";
          for(int c=0;c<A.Ncomp();c++)
            {
              for(int d=0;d<A.Mcomp();d++)
                {
                  s << A(i,j,c,d) << " ";
                }
            }
        }
    }
  return s;
}
}
