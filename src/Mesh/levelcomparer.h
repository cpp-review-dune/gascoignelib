#ifndef __levelcomparer_h
#define __levelcomparer_h

#include "compareclass.h"
#include  "gascoigne.h"

/*---------------------------------------------------*/

class LevelComparer
{
  const HierarchicalMesh&   Mesh;
  const Gascoigne::IntVector& v;

  public:
  
  LevelComparer(const HierarchicalMesh& HM, const Gascoigne::IntVector& vv) : 
    Mesh(HM), v(vv) {};

  int size() const { return v.size(); }
  int operator[](int i) const 
    {
      if (Mesh.dimension()==2)
	return Mesh.quad(v[i]).level();
      else
	return Mesh.hex(v[i]).level();
    }
};

/*---------------------------------------------------*/

#endif
