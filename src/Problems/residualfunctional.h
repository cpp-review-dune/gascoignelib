#ifndef  __ResidualFunctional_h
#define  __ResidualFunctional_h

#include  "functional.h"
#include  <string>
#include  <set>
#include  "dirichletdata.h"

/*-----------------------------------------*/


class ResidualFunctional : public Functional
{
protected:

  typedef  nvector<double>           Vector;

  int            _comp;
  std::set<int>  _col;
  double         _scale;

  DirichletData*   _DD;

public:

  ResidualFunctional();
  ResidualFunctional(const std::vector<std::string>& args);
  ~ResidualFunctional();

  void Construct(const std::vector<std::string>& args);

  std::string GetName() const {return "ResidualFunctional";}

  int           GetComp()   const {return _comp;}
  std::set<int> GetColors() const {return _col;}
  double        GetScale()  const { return _scale;}

  const DirichletData* GetDirichletData() const {return _DD;}
};


#endif
