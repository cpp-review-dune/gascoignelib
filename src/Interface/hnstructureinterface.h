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

namespace Gascoigne
{
class HNStructureInterface
{
public:

//
////  Con(De)structor 
//

  HNStructureInterface() {}
  virtual ~HNStructureInterface() {}

  virtual void ReInit(const MeshInterface* m)=0;
  virtual void MatrixDiag(int ncomp, MatrixInterface& A) const=0;
  virtual void SparseStructureDiag(SparseStructure* S) const=0;

  virtual void CondenseHanging(IntVector& indices) const=0;
  virtual void CondenseHanging(EntryMatrix& E, IntVector& indices) const=0;
  virtual void CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const=0;

  virtual void Average(GlobalVector& u) const=0;
  virtual void Distribute(GlobalVector& u) const=0;
  virtual void Zero(GlobalVector& u) const=0;
  virtual bool ZeroCheck(const GlobalVector& u) const=0;

  virtual int nhnodes() const =0;
};
}

#endif
