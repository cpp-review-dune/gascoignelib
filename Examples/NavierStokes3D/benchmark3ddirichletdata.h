#ifndef  __BenchMark3dDirichletData_h
#define  __BenchMark3dDirichletData_h

/////////////////////////////////////////////
///
///@brief
///  ... comments BenchMarkDirichletData
///
///
/////////////////////////////////////////////

#include  <string>
#include  "dirichletdata.h"

class BenchMark3dDirichletData : public DirichletData
{
protected:

  double vmax;

public:

//
///  Constructor 
//
  BenchMark3dDirichletData(const std::string& paramfile);
  std::string GetName() const {return "BenchMark3dDirichletData";}
  void operator()(Vector& b, const Vertex3d& v, int col) const;
};


#endif
