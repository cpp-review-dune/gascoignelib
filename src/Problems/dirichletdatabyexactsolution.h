#ifndef  __DirichletDataByExactSolution_h
#define  __DirichletDataByExactSolution_h

#include  "dirichletdata.h"
#include  "exactsolution.h"

/*-----------------------------------------*/


class DirichletDataByExactSolution : public DirichletData
{
protected:

  const ExactSolution* ES;

public:

  DirichletDataByExactSolution(const ExactSolution* es) 
    : ES(es), DirichletData() {assert(es);}

  std::string GetName() const {return "ExactSolution";}
  
  void operator()(Vector& b, const Vertex2d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }
  void operator()(Vector& b, const Vertex3d& v, int col)const{
    for(int c=0;c<b.size();c++) b[c] = (*ES)(c,v);
  }

};


#endif
