#include  "nodesparsestructureadaptor.h"

using namespace std;

/*-----------------------------------------*/

void NodeSparseStructureAdaptor::FillStencil(ColumnDiagStencil& S) const
{
  S.start(0) = 0;
  for(int i=0;i<SSP->n();i++)
    {
      int srowsize = _ncomp*SSP->rowsize(i);
      for(int c=0;c<_ncomp;c++)
	{
	  int ii = index(i,c);
	  int first = S.start(ii);
// 	  S.start(ii+1) = first + srowsize;
	  S.stop(ii) = first + srowsize;
	  
	  int id = 0;
	  for(SparseStructure::const_iterator p=SSP->rowbegin(i);
	      p!=SSP->rowend(i);p++)
	    {
	      for(int d=0;d<_ncomp;d++)
		{
		  S.col(first+id) = index(*p,d);
		  id++;
		}
	    }
	}
    }
}

/*-----------------------------------------*/

nvector<int> NodeSparseStructureAdaptor::GetIndicesDirichlet(int inode, const vector<int>& cv) const
{
  nvector<int> indices(cv.size());
  for(int ic=0;ic<cv.size();ic++)
    {
      indices[ic] = index(inode,cv[ic]);
    }
  return indices;
}
