#ifndef  __BasicIntegrator_h
#define  __BasicIntegrator_h



/////////////////////////////////////////////
////
////@brief
////  ... comments BasicIntegrator

////
////
/////////////////////////////////////////////


#include  "gascoigne.h"
#include  "integratorinterface.h"

class BasicIntegrator : public IntegratorInterface
{
 private:
  
  
 protected:
  
  mutable Gascoigne::FemFunction   NNN;
  mutable Gascoigne::TestFunction  NN, MM;
  mutable Gascoigne::FemFunction   UH;
  mutable Gascoigne::FemData       QH;
  
  void  universal_point(const FemInterface& FEM, Gascoigne::FemFunction& UH, const Gascoigne::LocalVector& U) const;
  void  universal_point(const FemInterface& FEM, Gascoigne::FemData& QH, const Gascoigne::LocalNodeData& Q) const;
  void  universal_point(Gascoigne::FemFunction& UH, const Gascoigne::LocalVector& U, const Gascoigne::FemFunction& NN) const;
  
 public:
  
  
  //
  ////  Con(De)structor 
  //

  BasicIntegrator();
  ~BasicIntegrator() {}
  
};


#endif
