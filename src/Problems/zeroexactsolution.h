#ifndef  __ZeroExactSolution_h
#define  __ZeroExactSolution_h

#include  "exactsolution.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments ZeroExactSolution

///
///
/////////////////////////////////////////////

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
}

#endif
