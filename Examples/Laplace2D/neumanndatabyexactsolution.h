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

class NeumannDataByExactSolution : public Gascoigne::NeumannData
{
private:

  const Gascoigne::Equation*      _EQ;
  const Gascoigne::ExactSolution* _ES;

protected:


public:


//
////  Con(De)structor 
//
  
  NeumannDataByExactSolution(const Gascoigne::Equation* eq, const Gascoigne::ExactSolution* es)
    : NeumannData(), _EQ(eq), _ES(es) { assert(es); assert(eq); }
  ~NeumannDataByExactSolution() {}

  std::string GetName() const {return "NeumannDataByExactSolution";}
  int GetNcomp() const { return _EQ->ncomp();}
  
  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Gascoigne::Vertex2d& v, const Gascoigne::Vertex2d& n, int col) const{
    b[0] += ( _ES->x(0,v)*n.x()+_ES->y(0,v)*n.y() ) * N.m();
  }
  
 };


#endif
