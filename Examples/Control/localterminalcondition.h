#ifndef  __LocalTerminalCondition_h
#define  __LocalTerminalCondition_h



/////////////////////////////////////////////
///
///@brief
///  ... comments LocalTerminalCondition

///
///
/////////////////////////////////////////////


#include  "domainrighthandside.h"
#include  "localequation.h"

class LocalTerminalCondition : public Gascoigne::DomainRightHandSide
{
public:


private:

  mutable const Gascoigne::FemFunction* q;

protected:

public:


//
///  Constructor 
//
  LocalTerminalCondition() {}
  
  std::string GetName() const { return "Local";}
  int GetNcomp() const {return 1;}  

  void SetFemData(Gascoigne::FemData& Q) const;
  double operator()(int c, const Gascoigne::Vertex2d& v) const;
};


#endif
