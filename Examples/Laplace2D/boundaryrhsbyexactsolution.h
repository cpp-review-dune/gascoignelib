#ifndef  __BoundaryRightHandSideByExactSolution_h
#define  __BoundaryRightHandSideByExactSolution_h


/////////////////////////////////////////////
////
////@brief
////  ... comments BoundaryRightHandSideByExactSolution

////
////
/////////////////////////////////////////////

#include  "boundaryrighthandside.h"
#include  "exactsolution.h"
#include  "equation.h"

class BoundaryRightHandSideByExactSolution : public Gascoigne::BoundaryRightHandSide
{
private:

  const Gascoigne::Equation*      _EQ;
  const Gascoigne::ExactSolution* _ES;

protected:


public:


//
////  Con(De)structor 
//
  
  BoundaryRightHandSideByExactSolution(const Gascoigne::Equation* eq, const Gascoigne::ExactSolution* es)
    : BoundaryRightHandSide(), _EQ(eq), _ES(es) { assert(es); assert(eq); }
  ~BoundaryRightHandSideByExactSolution() {}

  std::string GetName() const {return "BoundaryRightHandSideByExactSolution";}
  int GetNcomp() const { return _EQ->GetNcomp();}
  
  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Gascoigne::Vertex2d& v, const Gascoigne::Vertex2d& n, int col) const{
    b[0] += ( _ES->x(0,v)*n.x()+_ES->y(0,v)*n.y() ) * N.m();
  }
  
 };


#endif
