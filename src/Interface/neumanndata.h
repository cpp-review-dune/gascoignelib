#ifndef  __neumanndata_h
#define  __neumanndata_h

#include  "vertex.h"
#include  <set>
#include  "nvector.h"
#include  <string>
#include  "derivativevector.h"
#include  "gascoigne.h"
#include  "application.h"


//////////////////////////////////////////////
///
///@brief
/// Interface class for Neumann Boundary Conditions

/// void operator()(Vector& b, const Vertex2d& v, int col)
/// gets the coordinate v and color of boundarypart "col" and 
/// sets the values of b. b is a vector of length ncomp
///
//////////////////////////////////////////////

/*-----------------------------------------*/

class NeumannData : public Application
{
protected:

public:

  NeumannData() : Application() {}
  ~NeumannData() {}

  virtual std::string GetName() const=0;
  
  virtual void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex2d& v, const Vertex2d& n, int col) const {assert(0);}
  virtual void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex3d& v, const Vertex3d& n, int col) const {assert(0);}
};


#endif
