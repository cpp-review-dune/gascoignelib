#ifndef  __ZeroExactSolution_h
#define  __ZeroExactSolution_h



/////////////////////////////////////////////
///
///@brief
///  ... comments ZeroExactSolution

///
///
/////////////////////////////////////////////


#include  "exactsolution.h"


class ZeroExactSolution : public ExactSolution
{
public:


private:


protected:


public:


//
///  Constructor 
//

  ZeroExactSolution() {}
  std::string GetName() const{return "Zero";}
  double operator()(int c, const Vertex2d& v)const {return 0.;}

};


#endif
