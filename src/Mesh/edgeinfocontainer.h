#ifndef __EdgeInfoContainer_h
#define __EdgeInfoContainer_h

#include "edgeinfo.h"
#include "hierarchicalmesh2d.h"
#include "hierarchicalmesh3d.h"
#include "nvector.h"

/**********************************************************/

template<int DIM>
class EdgeInfoContainer : public nvector<EdgeInfo<DIM>*>
{

 protected:

  const HierarchicalMesh* _HMP;
  int                     _ncomp;

 public:

  EdgeInfoContainer<DIM>() {}
  ~EdgeInfoContainer<DIM>();

  void basicInit(const HierarchicalMesh*, int);
  void modifyHanging();

  const HierarchicalMesh* getMesh() const;

  void showStatistics() const;
};

#endif
