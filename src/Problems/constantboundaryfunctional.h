#ifndef  __ConstantBoundaryFunctional_h
#define  __ConstantBoundaryFunctional_h


#include  "boundaryfunctional.h"

/*-----------------------------------------*/


class ConstantBoundaryFunctional : public BoundaryFunctional
{
protected:

  int              comp;
  std::set<int>    colors;
  double           value;

public:


  ConstantBoundaryFunctional();
  ConstantBoundaryFunctional(const Equation& EQ, const std::vector<std::string>& args);
  ~ConstantBoundaryFunctional();
  void Construct(const Equation& EQ, const std::vector<std::string>& args);
  
  std::string GetName() const {return "ConstantBoundaryFunctional";}

  std::set<int> GetColors() const {return colors;}

  void AddColor(int    c) {colors.insert(c);}
  void SetComp (int    c) {comp =c;}
  void SetValue(double v) {value=v;}

  double J(const FemFunction& U, const Vertex2d& v) const;

};


#endif
