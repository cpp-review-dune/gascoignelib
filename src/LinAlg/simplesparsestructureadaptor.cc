#include  "simplesparsestructureadaptor.h"

/*-----------------------------------------*/

void SimpleSparseStructureAdaptor::FillStencil(ColumnDiagStencil& S) const
{
  S.start(0) = 0;
  for(int i=0;i<SSP->n();i++)
    {
      int srowsize = SSP->rowsize(i);
      int first = S.start(i);
      S.stop(i) = first + srowsize;
	  
      int id = 0;
      for(SparseStructure::const_iterator p=SSP->rowbegin(i);
	  p!=SSP->rowend(i);p++)
	{
	  S.col(first+id++) = *p;
	}
    }
}
