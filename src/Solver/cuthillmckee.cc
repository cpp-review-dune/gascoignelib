#include "cuthillmckee.h"
#include "ilupermutate.h"

using namespace std;

/* --------------------------------------------------------------- */

extern "C" void METIS_NodeND
(int*, int*, int*, int*, int*, int*, int*);

/* --------------------------------------------------------------- */

CuthillMcKee::CuthillMcKee (const StencilInterface *s)
{
  S=dynamic_cast<const ColumnStencil*>(s);
  assert(S);
//   M=m;
//   dimension=0;
}

/* --------------------------------------------------------------- */

CuthillMcKee::CuthillMcKee ()
{
//   M=0; 
  S=0;
//   dimension=0;
}

// /* --------------------------------------------------------------- */

// void CuthillMcKee::Permutate    (nvector<int> &perm, const Vertex2d v)
// {
//   dimension=2;
// //   dir2d=v;
//   Permutate(perm);
// }

// /* --------------------------------------------------------------- */

// void CuthillMcKee::Permutate    (nvector<int> &perm, const Vertex3d v)
// {
//   dimension=3;
// //   dir3d=v;
//   Permutate(perm);
// }

/* --------------------------------------------------------------- */

void CuthillMcKee::Permutate (nvector<int> &perm)
{
  // mit metis graph aufbauen
  //
  int n = S->n();

  perm.resize(n);

  vector<int> adj(n+1,0);
  vector<int> adjncy;
  
  int count=0;
  int c;
  adj[0]=0;
  for (int r=0;r<n;++r)
    {
      for (int p=S->start(r); p<S->stop(r);++p)
	{
	  c = S->col(p);
	  if (r==c) continue;
	  ++count;
	  adjncy.push_back(c);
	}
      adj[r+1]=count;
    }
  int numflag = 0;
  int options[8];
  options[0]=1;
  options[1]=3;
  options[2]=1;
  options[3]=1;
  options[4]=0;
  options[5]=1;
  options[6]=0;
  options[7]=1;
  vector<int> iperm(n);

  METIS_NodeND(&n,&adj[0],&adjncy[0],&numflag,&options[0],&perm[0],&iperm[0]);
  
//    assert((dimension==0)||(M->dimension()==dimension));
//      				   // Liste mit Nachbarn aufbauen
//    typedef SparseStructure::const_iterator SetIterator;
//    int n = S->n();
//    neighbors.resize(n);
//    for (unsigned int i=0;i<n;++i)
//      {
//        neighbors[i]=S->rowsize(i);
//      }
  
  
//  				   // Startvektor finden
//  //  stable_sort(perm.begin(),perm.end(),*this);
//    int starting_node=0;
//    for (int i=1;i<perm.size();++i)
//      if (operator()(perm[i],perm[starting_node])) starting_node=i;
  
//  				   // jetzt werden die Knoten so sortiert
//  				   // dass immer die Nachbarn eines Knotens
//  				   // kommen, sortiert nach neighbors
//    nvector<bool> versorgt(n,false);
//    int index = 0;
//    perm[index]=starting_node; versorgt[perm[index]]=true;
//    int st=0;
//    int en=0;
//    int st1=0;
//    while (en<n-1)
//      {
//        st1=index+1;
//        if (en<st) abort();
//        for (int i=st;i<=en;++i)
//  	{
//    	  for (int p=S->start(perm[i]);
//  	       p<S->stop(perm[i]);++p)
//    	    if (!versorgt[S->col(p)])
//    	      {
//    		++index; perm[index]=S->col(p);
//  		versorgt[perm[index]]=true;
//    	      }
//    	  if (index>st1)
//  	    stable_sort(perm.begin()+st1,perm.begin()+index+1,*this);
//  	}
//        st=en+1;
//        en=index;
//      }

//    vector<int> iperm(n);
//    for (int i=0;i<n;++i) iperm[perm[i]]=i;
  
//    char s[10];
//    sprintf(s,"l_%d",n);
//    ofstream aus(s);
//    for (int r=0;r<n;++r)
//      {
//        int row = perm[r];
//        for (int p = S->start(row);p!=S->stop(row);++p)
//  	aus << r << "\t" << iperm[S->col(p)] << endl;
//      }
  
//    aus.close();

}

// /* --------------------------------------------------------------- */

// bool CuthillMcKee::operator()(int i,int j) const
// {
//   if (neighbors[i]<neighbors[j]) return 1;
//   if (neighbors[i]>neighbors[j]) return 0;
//   if (dimension==2) return (((M->vertex2d(j)-M->vertex2d(i))*dir2d)>0);
//   if (dimension==3) return (((M->vertex3d(j)-M->vertex3d(i))*dir3d)>0);
//   return 0;
// }

