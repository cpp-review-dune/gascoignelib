#ifndef __ilupermutate_H
#define __ilupermutate_H

#include  "meshinterface.h"
#include  "compvector.h"
#include  "columnstencil.h"

class StreamDirection 
{
    int      dimension;
    int      dx,dy,dz;
    const MeshInterface* M;
    const ColumnStencil* S;
    const CompVector<double>&  X;

    void Permutate    (nvector<int> &perm);
    
    
  public:
    StreamDirection (const MeshInterface* m, const StencilInterface *s,
		     const CompVector<double>& x);
    
    void Permutate    (nvector<int> &perm,const nvector<int> d);
    
    bool operator() (int i,int j) const;
    double est      (int i,int j) const;
};

class VecDirection 
{
    Vertex2d dir2d;
    Vertex3d dir3d;
    int      dimension;
    const MeshInterface* M;

    void Permutate    (nvector<int> &perm);
    
  public:
    VecDirection (const MeshInterface* m);

    void Permutate    (nvector<int> &perm,nvector<double> v);

    bool operator()(int i,int j) const;
};

#endif








