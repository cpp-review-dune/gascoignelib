#include  "hnstructureq22d.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

HNStructureQ22d::HNStructureQ22d() : HNStructureQ12d()
{
  wei[0] =  0.375; 
  wei[1] =  0.75; 
  wei[2] = -0.125;

  lnoe[0][0]=0; lnoe[0][1]=2; lnoe[0][2]=1;
  lnoe[1][0]=0; lnoe[1][1]=6; lnoe[1][2]=3;
  lnoe[2][0]=2; lnoe[2][1]=8; lnoe[2][2]=5;
  lnoe[3][0]=6; lnoe[3][1]=8; lnoe[3][2]=7;
}

/*-----------------------------------------*/

void HNStructureQ22d::Average(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      const fixarray<3,int>& f = p->second;
      u.equ_node(p->first, wei[0], f[0], wei[1], f[1], wei[2], f[2]);
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::Distribute(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;
      for (int j=0; j<3; j++)
	{
	  u.add_node(p->second[j],wei[j],i);
	}
      u.zero_node(i);
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[2*ii+1];
      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices[p[0]]==f[1]) && (indices[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices[p[2]] = f[2];

      E.add_column     (p[0],p[2],weight(0));
      E.add_column     (p[1],p[2],weight(1));
      E.multiply_column(p[2],     weight(2));
      
      E.add_row        (p[0],p[2],weight(0));
      E.add_row        (p[1],p[2],weight(1));
      E.multiply_row   (p[2],     weight(2));
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHanging(IntVector& indices) const
{
  for(int ii=0; ii<4; ii++)
    {
      int j = lnoe[ii][2];
      int i = indices[j];

      if (hanging(i))
	indices[j] = regular_nodes(i)[2];
    }
}


}
