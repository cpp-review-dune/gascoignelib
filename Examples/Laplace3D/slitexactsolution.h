#ifndef  __SlitExactSolution_h
#define  __SlitExactSolution_h

#include  "exactsolution.h"
#include  "stringutil.h"
#include  "fadamath.h"


/////////////////////////////////////////////
///
///@brief
///  ... comments SlitExactSolution

///
///
/////////////////////////////////////////////




class SlitExactSolution : public ExactSolution
{
private:

public:

  std::string GetName() const {return "Slit";}


//
///  Constructor 
//
  SlitExactSolution() {}

  double operator()(int c, const Vertex2d& v)const 
  {
    double x = v.x();
    double y = v.y();
    double r = sqrt(x*x+y*y);

    double pi = GascoigneMath::pi();
    double theta;

    double fx = fabs(x);
    double fy = fabs(y);
    if(fx)
      {
	theta = atan(fy/fx);

	if     ( (x<0)&&(y>=0)) theta = pi-theta;
	else if( (x<0)&&(y<0))  theta += pi;
	else if( (x>0)&&(y<0))  theta = 2.*pi-theta;
      }
    else
      {
	if(y>=0) theta = 0.5*pi;
	else     theta = 1.5*pi;
      }
    return pow(r,0.5)*sin(0.5*theta);
  }
//   double x(int c, const Vertex2d& v)const
//   {
//     double x0 = v.x();
//     double y0 = v.y();
//     double r = sqrt(x0*x0+y0*y0)+1e-6;
//     return -y0/(r*r);
//   }
//   double y(int c, const Vertex2d& v)const
//   {
//     double x0 = v.x();
//     double y0 = v.y();
//     double r = sqrt(x0*x0+y0*y0)+1e-6;
//     return x0/(r*r);
//   }

};


#endif
