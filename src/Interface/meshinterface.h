#ifndef  __MeshInterface_h
#define  __MeshInterface_h


/////////////////////////////////////////////
///
///@brief
///  ... comments MeshInterface

///
///
/////////////////////////////////////////////

#include  "vertex.h"
#include  <set>
#include  <string>
#include  "gascoigne.h"

class MeshInterface
{
public:

  MeshInterface() {}
  virtual ~MeshInterface() {}

  virtual std::string GetName() const=0;

  virtual int  dimension() const=0;
  virtual int  nnodes()    const=0;
  virtual int  ncells()    const=0;

  virtual int  nodes_per_cell(int i)         const { assert(0);}
  virtual int  vertex_of_cell(int i, int ii)    const=0;
  virtual const Vertex2d& vertex2d(int i)  const { assert(0);}
  virtual const Vertex3d& vertex3d(int i)  const { assert(0);} 
  virtual std::set<int> GetColors()        const { assert(0);}
  virtual Gascoigne::IntVector  IndicesOfCell(int iq) const { assert(0);}
  virtual const Gascoigne::IntVector& Vertexo2n()     const { assert(0);}

  virtual const Gascoigne::IntVector& CellOnBoundary(int color)   const { assert(0);}
  virtual const Gascoigne::IntVector& LocalOnBoundary(int color)  const { assert(0);}

  virtual int VtkType(int i) const {assert(0);}

  // aus vtkCellTypes.h

// // Linear cells#define VTK_EMPTY_CELL     0
// #define VTK_VERTEX         1
// #define VTK_POLY_VERTEX    2
// #define VTK_LINE           3
// #define VTK_POLY_LINE      4
// #define VTK_TRIANGLE       5
// #define VTK_TRIANGLE_STRIP 6
// #define VTK_POLYGON        7
// #define VTK_PIXEL          8
// #define VTK_QUAD           9
// #define VTK_TETRA         10
// #define VTK_VOXEL         11
// #define VTK_HEXAHEDRON    12
// #define VTK_WEDGE         13
// #define VTK_PYRAMID       14

// // Quadratic, isoparametric cells
// #define VTK_QUADRATIC_EDGE       21
// #define VTK_QUADRATIC_TRIANGLE   22
// #define VTK_QUADRATIC_QUAD       23
// #define VTK_QUADRATIC_TETRA      24
// #define VTK_QUADRATIC_HEXAHEDRON 25

// // Special class of cells formed by convex group of points
// #define VTK_CONVEX_POINT_SET 41

// // Higher order cells in parametric form
// #define VTK_PARAMETRIC_CURVE        51
// #define VTK_PARAMETRIC_SURFACE      52
// #define VTK_PARAMETRIC_TRI_SURFACE  53
// #define VTK_PARAMETRIC_QUAD_SURFACE 54
// #define VTK_PARAMETRIC_TETRA_REGION 55
// #define VTK_PARAMETRIC_HEX_REGION   56}
};


#endif
