#ifndef  __ResidualFunctional_h
#define  __ResidualFunctional_h

#include  "functional.h"
#include  <string>
#include  <set>
#include  "dirichletdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class ResidualFunctional : public virtual Functional
{
protected:

  int            _comp;
  std::set<int>  _col;
  double         _scale;

  const DirichletData*   _DD;

public:

  ResidualFunctional();
  ~ResidualFunctional();
  ResidualFunctional(const ResidualFunctional& F) : Functional(F)
    {
       _comp  = F.GetComp();
       _col   = F.GetColors();
       _scale = F.GetScale();
       _DD    = F.GetDirichletData();
    }

  std::string GetName() const {return "ResidualFunctional";}

  int           GetComp()   const {return _comp;}
  int&          GetComp()         {return _comp;}
  
  std::set<int> GetColors() const {return _col;}
  std::set<int>&GetColors()       {return _col;}
  
  double        GetScale()  const { return _scale;}
  double&       GetScale()        { return _scale;}

  const DirichletData* GetDirichletData()   const {return _DD;}
  const DirichletData*& GetDirichletDataPointer() { return _DD;}
  
};
}

#endif
