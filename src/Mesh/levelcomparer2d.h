#ifndef __levelcomparer2d_h
#define __levelcomparer2d_h

#include "compareclass.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class LevelComparer2d
{
  const HierarchicalMesh2d&   Mesh;
  const IntVector& v;

  public:
  
  LevelComparer2d(const HierarchicalMesh2d& HM, const IntVector& vv) : 
    Mesh(HM), v(vv) {};

  int size() const { return v.size(); }
  int operator[](int i) const 
    {
      return Mesh.quad(v[i]).level();
    }
};
}

/*---------------------------------------------------*/

#endif
