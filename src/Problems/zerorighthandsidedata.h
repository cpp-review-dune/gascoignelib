#ifndef  __ZeroRightHandSideData_h
#define  __ZeroRightHandSideData_h

/////////////////////////////////////////////
///
///@brief
///  ... comments ZeroRightHandSideData

///
///
/////////////////////////////////////////////


#include  "righthandsidedata.h"


class ZeroRightHandSideData : public RightHandSideData
{
  int ncomp;

public:

//
///  Constructor 
//
  ZeroRightHandSideData(int nc) { ncomp = nc;}
  std::string GetName() const {return "zero";} 
  double operator()(int c, const Vertex2d& v)const {return 0.;}
  double operator()(int c, const Vertex3d& v)const {return 0.;}
  int GetNcomp() const { return ncomp;}
};


#endif
