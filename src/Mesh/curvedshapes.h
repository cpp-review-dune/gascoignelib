#ifndef __curvedshapes_h
#define __curvedshapes_h

#include <map>
#include <set>
#include <vector>
#include  <string>
#include "boundaryfunction.h"
#include "boundaryline.h"
#include "boundaryquad.h"
#include "vertex.h"

/******************************************************/

class CurvedShapes
{
 protected:
  
  std::set<int>                       colors;
  std::map<int,BoundaryFunction<2>* > ShapeOfColor2d;
  std::map<int,BoundaryFunction<3>* > ShapeOfColor3d;
  
 public:

  CurvedShapes();
  ~CurvedShapes();

  const BoundaryFunction<2>& GetShape2d(int col) const { return *(ShapeOfColor2d.find(col)->second);}
  const BoundaryFunction<3>& GetShape3d(int col) const { return *(ShapeOfColor3d.find(col)->second);}

  void newton(int col,Vertex2d& V) { ShapeOfColor2d[col]->newton(V); }
  void newton(int col,Vertex3d& V) { ShapeOfColor3d[col]->newton(V); };

  virtual void BasicInit(const std::vector<std::string>& names);
  void ReInit(const std::vector<BoundaryLine>& blines);
  void ReInit(const std::vector<BoundaryQuad>& bquads);

  int Curved(int col) const;
  int empty() const { return colors.size()==0;}
};

/******************************************************/

#endif
