#ifndef  __LocalMeshInterpretor_h
#define  __LocalMeshInterpretor_h



/////////////////////////////////////////////
////
////@brief
////  ... comments LocalMeshInterpretor

////
////
/////////////////////////////////////////////


#include  "q1gls2d.h"

class LocalMeshInterpretor : public Q1Gls2d
{
private:


protected:


public:


//
////  Con(De)structor 
//

  LocalMeshInterpretor();
  ~LocalMeshInterpretor() {}
  
  std::string GetName() const {return "Local";}

};


#endif
