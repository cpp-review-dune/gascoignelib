#include  "mginterpolatormatrix.h"


using namespace std;

/*-----------------------------------------*/
  
namespace Gascoigne
{
void MgInterpolatorMatrix::restrict_zero(GlobalVector& uL, const GlobalVector& ul) const
{
  uL.zero();
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  uL.add_node(ST.col(pos),val[pos],i,ul);
	}
    }
}

/*-----------------------------------------*/
  
void MgInterpolatorMatrix::prolongate_add(GlobalVector& ul, const GlobalVector& uL) const
{
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  ul.add_node(i,val[pos],ST.col(pos),uL);
	}
    }
}

/*-----------------------------------------*/

void MgInterpolatorMatrix::SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const 
{
  uL.zero();
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  if(val[pos]==1.)
	    uL.add_node(ST.col(pos),val[pos],i,ul);
	}
    }
}
}
