#ifndef __levelcomparer2d_h
#define __levelcomparer2d_h

#include "compareclass.h"

/*---------------------------------------------------*/

class LevelComparer2d
{
  const HierarchicalMesh2d&   Mesh;
  const nvector<int>& v;

  public:
  
  LevelComparer2d(const HierarchicalMesh2d& HM, const nvector<int>& vv) : 
    Mesh(HM), v(vv) {};

  int size() const { return v.size(); }
  int operator[](int i) const 
    {
      return Mesh.quad(v[i]).level();
    }
};

/*---------------------------------------------------*/

#endif
