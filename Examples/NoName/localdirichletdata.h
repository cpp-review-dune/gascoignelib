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
#include  "paramfile.h"

class LocalDirichletData : public DirichletData
{
protected:

  double vmax;

public:

//
///  Constructor 
//
  LocalDirichletData(const Gascoigne::ParamFile* paramfile);
  std::string GetName() const {return "Bench";}
  void operator()(Vector& b, const Vertex2d& v, int col) const;
};


#endif
