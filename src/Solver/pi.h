#ifndef  __Pi_h
#define  __Pi_h

#include  <map>
#include  "fixarray.h"
#include  "compvector.h"
#include  "gascoignemesh2d.h"
#include  "gascoignemesh3d.h"

namespace Gascoigne
{
/*-----------------------------------------*/

class Pi
{
 protected:
  
  std::map<int,fixarray<2,int> > edge;
  std::map<int,fixarray<4,int> > face;
  std::map<int,fixarray<8,int> > cell;

  void Init2d(const GascoigneMesh2d* MP);
  void Init3d(const GascoigneMesh3d* MP);
  
 public:
  
  Pi();

  void Init(const MeshInterface* MP);

  void vmult(CompVector<double>& y, const CompVector<double>& x, 
	     double s=1.) const;
};
}

#endif
