#ifndef  __DirichletDataByColor_h
#define  __DirichletDataByColor_h

#include  "dirichletdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
  class DirichletDataByColor : public DirichletData
    {
    protected:

      std::set<int>   __cols;
      nvector<int>    __comps;
      nvector<double> __scales;

    public:
      
      DirichletDataByColor(nvector<int> comps,const std::set<int>& cl, nvector<double> s);
      DirichletDataByColor(int comps, std::set<int>& cl, double s);
      DirichletDataByColor(const std::vector<std::string>& args);
      ~DirichletDataByColor() {}
      
      std::string GetName() const {return "DirichletDataByColor";}

      std::set<int> preferred_colors()const;
      void operator()(DoubleVector& b, const Vertex2d& V,int color) const;
      void operator()(DoubleVector& b, const Vertex3d& V,int color) const;
    };
}

#endif
