#include "deletecells.h"
#include "quad.h"
#include "hex.h"
#include "boundarycell.h"

/*---------------------------------------------------*/

template<class C>
void delete_cells(const IntSet& coarselist, 
		  std::vector<C>& liste, 
		  const nvector<int>& co2n, const nvector<int>& vo2n)
{
  for(unsigned oi=0;oi<co2n.size();++oi)
    {
      int ni = co2n[oi];
      if(ni>=0)  
	{
	  C q(liste[oi]);

	  for(unsigned i=0;i<q.vertex().size();++i)
	    {
	      q.vertex(i) = vo2n[q.vertex(i)];
	      if(q.vertex(i)==-1)
		{
		  std::cerr << "Vertex invalid in "<<oi<<" "<<ni<<std::endl;
		  //std::cerr << vo2n[liste[oi].vertex(i)]<<std::endl;
		  //std::cerr << q.vertex();
		  abort();
		}
	    }
	  IntSet::iterator  p=coarselist.find(oi);
	  if(p!=coarselist.end())
	    {
	      q.childs().resize(0);
	    }
	  int qf = -1;
	  if(q.father()!=-1) qf = co2n[q.father()]; 
	  q.father()= qf;
	  if(q.sleep())
	    {
	      for(int i=0;i<q.nchilds();++i)
		{
		  q.child(i) = co2n[q.child(i)];
		}
	    }
	  liste[ni] = q;
	}
    }
}

/*---------------------------------------------------*/

template void delete_cells<Quad>(const IntSet&, std::vector<Quad>&, 
				 const nvector<int>&, const nvector<int>&);
template void delete_cells<Hex>(const IntSet&, std::vector<Hex>&, 
				 const nvector<int>&, const nvector<int>&);

template void delete_cells<BoundaryCell<2> >(const IntSet&, std::vector<BoundaryCell<2> >&, 
				const nvector<int>&, const nvector<int>&);

template void delete_cells<BoundaryCell<4> >(const IntSet&, std::vector<BoundaryCell<4> >&, 
				const nvector<int>&, const nvector<int>&);
				
