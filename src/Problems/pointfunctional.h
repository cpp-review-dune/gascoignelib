#ifndef  __PointFunctional_h
#define  __PointFunctional_h

#include  <vector>
#include  "vertex.h"
#include  "functional.h"

/*-----------------------------------------*/


class PointFunctional : public Functional
{
protected:

  std::vector<Vertex2d>  v;
  nvector<double>        w;
  std::vector<int>       ids;
  std::string            type;

  int mycomp;
  
public:

  PointFunctional();
  PointFunctional(const Equation& EQ, const std::vector<std::string>& args);
  ~PointFunctional();

  void Init(const Equation& EQ, const std::vector<std::string>& args);
  
  std::string GetName() const {return "PointFunctional";}

  int size() const {return w.size();}
  const std::vector<Vertex2d>& points()      const { return v;}
  const nvector<double>&  weights()     const { return w;}
  const Vertex2d& point (int i)              const { return v[i];}
        Vertex2d& point (int i)                    { return v[i];}

  int GetComp() const { return mycomp; }
  
  double          weight(int i) const {return w[i];}
  int             id    (int i) const {return ids[i];}
};


#endif
