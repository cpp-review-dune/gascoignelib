#ifndef  __righthandsidedata_h
#define  __righthandsidedata_h

#include  "gascoigne.h"
#include  "vertex.h"
#include  <string>
#include  "nvector.h"
#include  "application.h"

//////////////////////////////////////////////
///
///@brief
/// Interface class for RightHandSideData

///
///
//////////////////////////////////////////////

/*-----------------------------------------*/

class RightHandSideData : public Application
{
protected:

public:

  RightHandSideData() : Application() {}
  ~RightHandSideData() {}

  virtual double operator()(int c, const Vertex2d& v) const {assert(0); return 0;}
  virtual double operator()(int c, const Vertex3d& v) const {assert(0); return 0;}

  virtual void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex2d& v) const 
    {
      for(int c=0;c<GetNcomp();c++)
	{
	  b[c] += N.m()* (*this)(c,v);
	}
    }
  virtual void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, const Vertex3d& v) const 
    {
      for(int c=0;c<GetNcomp();c++)
	{
	  b[c] += N.m()* (*this)(c,v);
	}
    }
/*   virtual const nvector<double>& GetWeights()   const { assert(0);  return ;} */
/*   virtual const std::vector<Vertex2d>& GetPoints2d() const { assert(0); } */
/*   virtual const std::vector<Vertex3d>& GetPoints3d() const { assert(0); } */
};

#endif
