#ifndef  __LocalDirichletData_h
#define  __LocalDirichletData_h

/////////////////////////////////////////////
///
///@brief
///  ... comments LocalDirichletData
///
///
/////////////////////////////////////////////

#include  <string>
#include  "dirichletdata.h"

class LocalDirichletData : public DirichletData
{
protected:

  double vmax;

public:

//
///  Constructor 
//
  LocalDirichletData(const std::string& paramfile);
  std::string GetName() const {return "Bench";}
  void operator()(Vector& b, const Vertex2d& v, int col) const;
};


#endif
