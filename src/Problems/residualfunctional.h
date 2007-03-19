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

      nvector<int>    __comps;
      nvector<double> __scales;
      std::set<int>   __cols;

      const DirichletData*   __DD;
  
    public:
  
      ResidualFunctional();
      ~ResidualFunctional();
      ResidualFunctional(const ResidualFunctional& F) : Functional(F)
	{
	  __comps  = F.GetComps();
	  __scales = F.GetScales();
	  __cols   = F.GetColors();
	  __DD     = F.GetDirichletData();
	}

      std::string GetName() const {return "ResidualFunctional";}

      nvector<int>     GetComps()  const {return __comps;}
      nvector<int>&    GetComps()        {return __comps;}

      nvector<double>  GetScales() const { return __scales;}
      nvector<double>& GetScales()       { return __scales;}
  
      std::set<int>    GetColors() const {return __cols;}
      std::set<int>&   GetColors()       {return __cols;}
  

      const DirichletData* GetDirichletData()   const { return __DD;}
      const DirichletData*& GetDirichletDataPointer() { return __DD;}
  
    };
}

#endif
