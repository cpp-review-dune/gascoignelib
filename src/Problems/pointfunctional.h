#ifndef  __PointFunctional_h
#define  __PointFunctional_h

#include  <vector>
#include  "vertex.h"
#include  "functional.h"

/*-----------------------------------------*/

class PointFunctional : public Functional
{
protected:

  std::vector<Vertex2d>  v2d;
  std::vector<Vertex3d>  v3d;
  nvector<double>        w;
  std::vector<int>       ids;
  std::string            type;

  int mycomp;
  
public:

  PointFunctional(const std::vector<std::string>& args);
  ~PointFunctional();

  void Init(const std::vector<std::string>& args);
  
  std::string GetName() const {return "PointFunctional";}

  int size() const {return w.size();}
  const std::vector<Vertex2d>& points2d()    const { return v2d;}
  const std::vector<Vertex3d>& points3d()    const { return v3d;}
  const nvector<double>&  weights()          const { return w;}

  int GetComp() const { return mycomp; }
};


#endif
