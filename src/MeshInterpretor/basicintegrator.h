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
#include  "derivativevector.h"

using namespace Gascoigne;

class BasicIntegrator : public IntegratorInterface
{
 private:
  
  
 protected:
  
  mutable FemFunction   NNN;
  mutable DerivativeVector  NN, MM;
  mutable FemFunction      UH;
  mutable FemData  QH;
  
  void  universal_point(const FemInterface& FEM, FemFunction& UH, const LocalVector& U) const;
  void  universal_point(const FemInterface& FEM, FemData& QH, const LocalData& Q) const;
  void universal_point(FemFunction& UH, const LocalVector& U, const FemFunction& NN) const;
  
 public:
  
  
  //
  ////  Con(De)structor 
  //

  BasicIntegrator();
  ~BasicIntegrator() {}
  
};


#endif
