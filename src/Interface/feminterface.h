#ifndef  __FemInterface_h
#define  __FemInterface_h

#include  "vertex.h"
#include  "nmatrix.h"
#include  "gascoigne.h"
#include  <string>

/*-----------------------------------------*/

class FemInterface
{
public:

  typedef nmatrix<double>    Matrix;

protected:

  void error(const std::string& s) const
    {
      std::cout << "FemInterface::" << s <<" not written\n"; 
      abort();
    }

public:

  FemInterface() {}
  virtual ~FemInterface() {}

  virtual std::string GetName() const=0;

  virtual int    n()         const=0;
/*   virtual double N   (int i) const=0; */
/*   virtual double N_x (int i) const=0; */
/*   virtual double N_y (int i) const=0; */
/*   virtual double N_z (int i) const { return 0.;} */
/*   virtual double N_D (int i) const { return 0.;} */
  virtual double J()         const=0;
  virtual double G()         const=0;

  virtual void x(Vertex2d& v) const {assert(0);}
  virtual void x(Vertex3d& v) const {assert(0);}

  virtual void normal(Vertex2d& v) const {assert(0);}
  virtual void normal(Vertex3d& v) const {assert(0);}

  virtual void point(const Vertex2d& v) const
    { error("point 2d"); }
  virtual void point(const Vertex3d& v) const
    { error("point 3d"); }

  virtual void  point_boundary(int ie, const Vertex1d& v) const
    { error("point_boundary 1d"); };
  virtual void  point_boundary(int ie, const Vertex2d& v) const
    { error("point_boundary 2d"); };

/*   virtual void init(const Matrix& M)=0; */
  virtual void ReInit(const Matrix& M) const=0;

  virtual void  init_test_functions(Gascoigne::TestFunction& Phi, double w, int i) const=0;
};


#endif
