#ifndef  __LocalDirichletData_h
#define  __LocalDirichletData_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalDirichletData

////
////
/////////////////////////////////////////////

#include  "dirichletdata.h"

class LocalDirichletData : public DirichletData
{
private:


protected:


public:


//
////  Con(De)structor 
//
  LocalDirichletData() {}
  ~LocalDirichletData() {}

  std::string GetName() const {return "Local";}

  void operator()(Vector& b, const Vertex2d& v, int col) const;

};


#endif
