#ifndef  __dirichletdata_h
#define  __dirichletdata_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>

//////////////////////////////////////////////
///
///@brief
/// Interface class for Dirichlet Boundary Conditions

/// void operator()(Vector& b, const Vertex2d& v, int col)
/// gets the coordinate v and color of boundarypart "col" and 
/// sets the values of b. b is a vector of length ncomp
///
//////////////////////////////////////////////


/*-----------------------------------------*/

class DirichletData
{
protected:

  typedef  nvector<double>         Vector;

public:

  DirichletData() {}
  virtual ~DirichletData() {}

  virtual std::string GetName() const=0;

  virtual void operator()(Vector& b, const Vertex2d& v, int col) const {}
  virtual void operator()(Vector& b, const Vertex3d& v, int col) const {}

  virtual std::set<int> preferred_colors()const {return std::set<int>();}

  virtual void SetTime(double time, double dt) {}
};


#endif
