#ifndef  __NeumannDataByExactSolution_h
#define  __NeumannDataByExactSolution_h


/////////////////////////////////////////////
////
////@brief
////  ... comments NeumannDataByExactSolution

////
////
/////////////////////////////////////////////

#include  "neumanndata.h"
#include  "exactsolution.h"
#include  "equation.h"

class NeumannDataByExactSolution : public NeumannData
{
private:

  const Equation*      _EQ;
  const ExactSolution* _ES;

protected:


public:


//
////  Con(De)structor 
//
  
  NeumannDataByExactSolution(const Equation* eq, const ExactSolution* es)
    : NeumannData(), _EQ(eq), _ES(es) { assert(es); assert(eq); }
  ~NeumannDataByExactSolution() {}

  std::string GetName() const {return "NeumannDataByExactSolution";}
  int GetNcomp() const { return _EQ->ncomp();}
  
  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex2d& v, const Vertex2d& n, int col) const{
    b[0] += ( _ES->x(0,v)*n.x()+_ES->y(0,v)*n.y() ) * N.m();
  }
  
 };


#endif
