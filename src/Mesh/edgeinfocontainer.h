#ifndef __EdgeInfoContainer_h
#define __EdgeInfoContainer_h

#include "edgeinfo.h"
#include "hierarchicalmesh2d.h"
#include "hierarchicalmesh3d.h"
#include "nvector.h"

/**********************************************************/

namespace Gascoigne
{

class EdgeInfoContainerInterface
{
 public:

  virtual const HierarchicalMesh* GetMesh() const=0;
  virtual void BasicInit(const HierarchicalMesh*, int)=0;
};

template<int DIM>
  class EdgeInfoContainer : public virtual EdgeInfoContainerInterface, public nvector<EdgeInfo<DIM>*>
{

 protected:

  const HierarchicalMesh* _HMP;
  int                     _ncomp;

 public:

  EdgeInfoContainer<DIM>() {}
  ~EdgeInfoContainer<DIM>();

  void BasicInit(const HierarchicalMesh*, int);
  void ModifyHanging();

  const HierarchicalMesh* GetMesh() const { return _HMP; }

  void ShowStatistics() const;
};
}

#endif
