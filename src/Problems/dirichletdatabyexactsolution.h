#ifndef  __DirichletDataByExactSolution_h
#define  __DirichletDataByExactSolution_h

#include  "dirichletdata.h"
#include  "exactsolution.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class DirichletDataByExactSolution : public DirichletData
{
protected:

  const ExactSolution* ES;

public:

  DirichletDataByExactSolution(const ExactSolution* es) 
    : ES(es), DirichletData() {assert(es);}

  std::string GetName() const {return "ExactSolution";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }
  void operator()(DoubleVector& b, const Vertex3d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }

};
}

#endif
