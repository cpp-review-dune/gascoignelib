#include "ilupermutate.h"

/* =============================================================== */
/* ============== VecDirection =================================== */
/* =============================================================== */

VecDirection::VecDirection (const MeshInterface* m)
{
  M=m;
  dimension=0;
}

/* --------------------------------------------------------------- */

void VecDirection::Permutate (nvector<int> &perm)
{
  assert(dimension==M->dimension());
  stable_sort(perm.begin(),perm.end(), *this);
}

/* --------------------------------------------------------------- */

void VecDirection::Permutate (nvector<int> &perm,nvector<double> v)
{
  if (v.size()==0) 
    {
      v.resize(M->dimension(),0);
      v[0]=1.;
    }
  
  dimension=v.size();
  if (dimension==2)
    {
      dir2d[0]=v[0];
      dir2d[1]=v[1];
    }
  else if (dimension==3)
    {
      dir3d[0]=v[0];
      dir3d[1]=v[1];
      dir3d[2]=v[2];
    }
  else assert(0);
  
  Permutate(perm);
}


/* --------------------------------------------------------------- */

bool VecDirection::operator()(int i,int j) const
{
  if (dimension==2)
    return (((M->vertex2d(j)-M->vertex2d(i))*dir2d)>0);
  if (dimension==3)
    return (((M->vertex3d(j)-M->vertex3d(i))*dir3d)>0);
  assert(0);
}

/* =============================================================== */
/* ============== StreamDirection ================================ */
/* =============================================================== */


StreamDirection::StreamDirection (const MeshInterface* m,
				  const StencilInterface *s,
				  const CompVector<double>& x):
  X(x)
{
  M=m;
  S=dynamic_cast<const ColumnStencil*>(s);
  assert(S);
  dimension=0;
  dx=dy=dz=0;
}

/* --------------------------------------------------------------- */

void StreamDirection::Permutate    (nvector<int> &perm)
{
  assert(dimension==M->dimension());
  assert(perm.size()==M->nnodes());
  assert(perm.size()==X.n());
  assert(X.ncomp()>dx);
  assert(X.ncomp()>dy);
  assert(X.ncomp()>dz);
  int n=X.n();
  Vertex2d h2;
  Vertex3d h3;
				   // Zu allen Knoten den Nachfolger
				   // und Vorgaenger finden.
  nvector<int> next(n,-1);
  nvector<int> prev(n,-1);
				   // Matrixzeilen durchlaufen
  for (int row=0;row<n;++row)
    {
      assert (next[row]==-1);
      int best = -1;
      double bv = -1;
      				       // spalten
      for (int p=S->start(row);p<S->stop(row);++p)
	{
					   // gleicher Knoten
	  int col = S->col(p);
	  if (col==row) continue;
					   // hat noch keinen Vorgaenger
	  if (prev[col]==-1)
					     // ist ganz prima
	    if (est(row,col)>bv)
	      {
		bv=est(row,col);
		best=col;
	      }
	}
      if (bv>0.9)
	{
	  next[row]=best;
	  prev[best]=row;
	}
    }

				   // Kreise Finden und tot machen
  for (int i=0;i<n;++i)
    {
      int k=i;
      while ((next[k]!=-1)&&(next[k]!=i)) k=next[k];
      if (next[k]==i)
	{
	  prev[next[k]]=-1;
	  next[k]=-1;
	}
    }
  
  
				   // erster Versuch:
				   //
				   // Sequenzen werden aneinandergehaengt
				   // es wird nicht darauf geachtet, in
				   // welcher Reihenfolge die Sequenzen kommen.
  
  
  
  int count=0;
  int index=0;
  std::set<int> start;
  for (int i=0;i<n;++i) if (prev[i]==-1)
    start.insert(i);
    
  
  for (std::set<int>::const_iterator i=start.begin();i!=start.end();++i)
    {
      int k=*i;
      do
	{
	  perm[index]=k;
	  ++index;
	  k=next[k];
	}
      while (k!=-1);
      
    }
}

/* --------------------------------------------------------------- */

void StreamDirection::Permutate    (nvector<int> &perm,const nvector<int> d)
{
  dimension=d.size();
  assert(dimension>1);
  dx=d[0];
  dy=d[1];
  if (dimension==3)
    dz=d[2];
  Permutate(perm);
}

/* --------------------------------------------------------------- */


bool StreamDirection::operator()(int i,int j) const
{
  if (dimension==2)
    {
      numfixarray<2,double> a = M->vertex2d(j)-M->vertex2d(i);
      double sc = a[0]*X(i,dx)+a[1]*X(i,dy);
      double n  = sqrt(a[0]*a[0]+a[1]*a[1])*sqrt(X(i,dx)*X(i,dx)+X(i,dy)*X(i,dy));
      return (sc>0.75*n);
    }
  if (dimension==3)
    {
      numfixarray<3,double> a = M->vertex3d(j)-M->vertex3d(i);
      double sc = a[0]*X(i,dx)+a[1]*X(i,dy)+a[2]*X(i,dz);
      double n  = sqrt(a[0]*a[0]+a[1]*a[1])*
		  sqrt(X(i,dx)*X(i,dx)+X(i,dy)*X(i,dy)+X(i,dz)*X(i,dz));
      return (sc>0.5*n);
    }
}

/* --------------------------------------------------------------- */


double StreamDirection::est(int i,int j) const
{
  if (dimension==2)
    {
      numfixarray<2,double> a = M->vertex2d(j)-M->vertex2d(i);
      double sc = a[0]*(X(i,dx)+X(j,dx))+a[1]*(X(i,dy)+X(j,dy));
      double n  = sqrt(a[0]*a[0]+a[1]*a[1])*
		  sqrt((X(i,dx)+X(j,dx))*(X(i,dx)+X(j,dx))+
		       (X(i,dy)+X(j,dy))*(X(i,dy)+X(j,dy)));
      if (n<0.000001) return 0.51;
      return sc/n;
    }
  if (dimension==3)
    {
      numfixarray<3,double> a = M->vertex3d(j)-M->vertex3d(i);
      double sc = a[0]*X(i,dx)+a[1]*X(i,dy)+a[2]*X(i,dz);
      double n  = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]+a[3])*
		  sqrt(X(i,dx)*X(i,dx)+X(i,dy)*X(i,dy)+X(i,dz)*X(i,dz));
            if (n<0.0000001) return 0.51;
      return sc/n;
    }
}

