#ifndef __ilupermutate_H
#define __ilupermutate_H

#include  "meshinterface.h"
#include  "compvector.h"
#include  "columnstencil.h"

namespace Gascoigne
{
class StreamDirection 
{
    int      dimension;
    int      dx,dy,dz;
    const MeshInterface* M;
    const ColumnStencil* S;
    const GlobalVector&  X;

    void Permutate    (IntVector &perm);
    
    
  public:
    StreamDirection (const MeshInterface* m, const StencilInterface *s,
		     const GlobalVector& x);
    
    void Permutate    (IntVector &perm,const IntVector d);
    
    bool operator() (int i,int j) const;
    double est      (int i,int j) const;
};

class VecDirection 
{
    Vertex2d dir2d;
    Vertex3d dir3d;
    int      dimension;
    const MeshInterface* M;

    void Permutate    (IntVector &perm);
    
  public:
    VecDirection (const MeshInterface* m);

    void Permutate    (IntVector &perm,DoubleVector v);

    bool operator()(int i,int j) const;
};
}

#endif








