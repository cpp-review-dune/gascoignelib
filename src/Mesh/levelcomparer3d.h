#ifndef __levelcomparer3d_h
#define __levelcomparer3d_h

#include  "compareclass.h"
#include  "gascoigne.h"

/*---------------------------------------------------*/

class LevelComparer3d
{
  const HierarchicalMesh3d&   Mesh;
  const Gascoigne::IntVector& v;

  public:
  
  LevelComparer3d(const HierarchicalMesh3d& HM, const Gascoigne::IntVector& vv) : 
    Mesh(HM), v(vv) {};

  int size() const { return v.size(); }
  int operator[](int i) const 
    {
      return Mesh.hex(v[i]).level();
    }
};

/*---------------------------------------------------*/

#endif
