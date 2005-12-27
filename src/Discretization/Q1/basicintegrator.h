#ifndef  __BasicIntegrator_h
#define  __BasicIntegrator_h


#include  "gascoigne.h"
#include  "integratorinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments BasicIntegrator

////
////
/////////////////////////////////////////////

class BasicIntegrator : public IntegratorInterface
{
 private:
  
  
 protected:
  
  mutable FemFunction   _NNN;
  mutable TestFunction  _NN;
  mutable FemFunction   _UH;
  mutable FemData       _QH;
  mutable CellData      _QCH;
  
  void  universal_point(const FemInterface& FEM, FemFunction& UH, const LocalVector& U) const;
  void  universal_point(CellFunction& UCH, const LocalCellVector& UC) const;
  void  universal_point(FemFunction& UH, const LocalVector& U, const FemFunction& NN) const;
  
  void  universal_point(const FemInterface& FEM, FemData& QH, const LocalNodeData& Q) const;
  void  universal_point(CellData& QCH, const LocalCellData& QC) const;

 public:
  
  
  //
  ////  Con(De)structor 
  //

  BasicIntegrator();
  ~BasicIntegrator() {}
  
};
}

#endif
