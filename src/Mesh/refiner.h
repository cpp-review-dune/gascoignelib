#ifndef __refiner_h
#define __refiner_h

#include "meshinterface.h"
#include "nmatrix.h"

/*************************************************************/

class PointRefiner
{
  const MeshInterface& H;
  const Vertex2d& V;

  nmatrix<double> A;
  bool VertexInQuad(int);

 public:

  PointRefiner(const MeshInterface& h, const Vertex2d& v) 
    : H(h), V(v), A(2,2) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class CircleRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;

  bool QuadOnRadius(int) const;

 public:

  CircleRefiner(const MeshInterface& h, const Vertex3d& v, double r) 
    : H(h), V(v), R(r) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class CylinderRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;
  int D;

  bool QuadInCylinder(int) const;

 public:

  CylinderRefiner(const MeshInterface& h, const Vertex3d& v, 
		  double r, int d)
    : H(h), V(v), R(r), D(d) {}

  void BuildCellList(std::vector<int>&);
};

/*************************************************************/

class BallRefiner 
{
  const MeshInterface&  H;
  const Vertex3d&  V;
  double R;
  
  bool QuadInBall(int) const;

 public:
  
  BallRefiner(const MeshInterface& h, const Vertex3d& v, 
	      double r)
    : H(h), V(v), R(r) {}

  void BuildCellList(std::vector<int>&);
};


/*************************************************************/

#endif
