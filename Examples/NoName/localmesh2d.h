#ifndef  __LocalMesh2d_h
#define  __LocalMesh2d_h

#include  "gascoignemesh2d.h"
#include  "boundaryfunction.h"
#include  "levelmesh2d.h"

#include  <map>

/*-----------------------------------------*/

class LocalMesh2d : public GascoigneMesh2d
{
protected:

  int _count;
  std::vector<nvector<double> > _q;

  std::map<int,BoundaryFunction<2>*>  MyShapes;

  void boundary_newton2d(){}

public:

  LocalMesh2d();
  void SetCoordinates(const LevelMesh2d* LM);

};

#endif






