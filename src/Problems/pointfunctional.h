#ifndef  __PointFunctional_h
#define  __PointFunctional_h

#include  <vector>
#include  "vertex.h"
#include  "functional.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class PointFunctional : public Functional
{
protected:

  std::vector<Vertex2d>  v2d;
  std::vector<Vertex3d>  v3d;
  DoubleVector        w;
  std::vector<int>       ids;
  std::string            type;

  int mycomp;
  
public:

  PointFunctional(const std::vector<std::string>& args);
  PointFunctional() {};
  PointFunctional(const PointFunctional& F) : Functional(F) 
    {
      v2d = F.points2d();
      v3d = F.points3d();
      w   = F.weights() ;
      ids = F.GetIds();
      type = F.GetType();

      mycomp = F.GetComp();
      
    };
  ~PointFunctional();

  void Init(const std::vector<std::string>& args);
  
  std::string GetName() const {return "PointFunctional";}

  int size() const {return w.size();}
  const std::vector<Vertex2d>& points2d()    const { return v2d;}
  const std::vector<Vertex3d>& points3d()    const { return v3d;}
  const DoubleVector&  weights()          const { return w;}
  const std::string& GetType()               const { return type;}
  const std::vector<int>& GetIds()           const { return ids;}

  int GetComp() const { return mycomp; }
};
}

#endif
