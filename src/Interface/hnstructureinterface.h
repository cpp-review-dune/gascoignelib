#ifndef  __HNStructureInterface_h
#define  __HNStructureInterface_h


/////////////////////////////////////////////
////
////@brief
////  ... comments HNStructureInterface

////
////
/////////////////////////////////////////////

#include  "gascoigne.h"
#include  "meshinterface.h"
#include  "matrixinterface.h"
#include  "sparsestructure.h"

/*-----------------------------------------*/

class HNStructureInterface
{
private:


protected:


public:


//
////  Con(De)structor 
//

  HNStructureInterface() {}
  virtual ~HNStructureInterface() {}

  virtual void ReInit(const MeshInterface* m)=0;
  virtual void MatrixDiag(int ncomp, MatrixInterface& A) const=0;
  virtual void SparseStructureDiag(SparseStructure* S) const=0;

  virtual void CondenseHanging(nvector<int>& indices) const=0;
  virtual void CondenseHanging(EntryMatrix& E, nvector<int>& indices) const=0;
  virtual void CondenseHangingPatch(EntryMatrix& E, nvector<int>& indices) const=0;

  virtual void Average(Gascoigne::GlobalVector& u) const=0;
  virtual void Distribute(Gascoigne::GlobalVector& u) const=0;
  virtual void Zero(Gascoigne::GlobalVector& u) const=0;
  virtual bool ZeroCheck(const Gascoigne::GlobalVector& u) const=0;

};


#endif
