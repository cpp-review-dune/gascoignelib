#ifndef  __ExactSolution_h
#define  __ExactSolution_h

//////////////////////////////////////////////
///
///@brief
/// Interface class for ExactSolution

///
///
//////////////////////////////////////////////


#include  "vertex.h"
#include  <string>
#include  "application.h"

/*-----------------------------------------*/

class ExactSolution : public Application
{
protected:

  double eps;

public:

  ExactSolution(): Application(), eps(1.e-6) {}

  virtual ~ExactSolution() {}

  virtual std::string GetName() const=0;

  ////////// 2d

  virtual double operator()(int c, const Vertex2d& v)const {return 0.;}

  virtual double x         (int c, const Vertex2d& v)const;
  virtual double y         (int c, const Vertex2d& v)const;
  virtual double z         (int c, const Vertex2d& v)const {assert(0);}
  virtual double xx        (int c, const Vertex2d& v)const;
  virtual double yx        (int c, const Vertex2d& v)const;
  virtual double xy        (int c, const Vertex2d& v)const;
  virtual double yy        (int c, const Vertex2d& v)const;

  ////////// 3d

  virtual double operator()(int c, const Vertex3d& v)const {return 0.;}

  virtual double x         (int c, const Vertex3d& v)const;
  virtual double y         (int c, const Vertex3d& v)const;
  virtual double z         (int c, const Vertex3d& v)const;
  virtual double xx        (int c, const Vertex3d& v)const;
  virtual double yy        (int c, const Vertex3d& v)const;
  virtual double zz        (int c, const Vertex3d& v)const;
  virtual double xy        (int c, const Vertex3d& v)const;
  virtual double yz        (int c, const Vertex3d& v)const;
  virtual double xz        (int c, const Vertex3d& v)const;
};


#endif
