#ifndef  __BenchMarkDirichletData_h
#define  __BenchMarkDirichletData_h

/////////////////////////////////////////////
///
///@brief
///  ... comments BenchMarkDirichletData
///
///
/////////////////////////////////////////////

#include  <string>
#include  "dirichletdata.h"

class BenchMarkDirichletData : public DirichletData
{
protected:

  double vmax;

public:

//
///  Constructor 
//
  BenchMarkDirichletData(const std::string& paramfile);
  std::string GetName() const {return "Bench";}
  void operator()(Vector& b, const Vertex2d& v, int col) const;
};


#endif
