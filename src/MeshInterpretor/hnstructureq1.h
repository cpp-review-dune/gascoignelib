#ifndef  __HNStructureQ1_h
#define  __HNStructureQ1_h

#include  "matrixinterface.h"
#include  "compvector.h"
#include  "meshinterface.h"
#include  "sparsestructure.h"
#include  <map>

/*-----------------------------------------*/

class HNStructureQ1
{
protected:

  typedef  CompVector<double>     Vector;
  typedef  fixarray<3,int>        EdgeVector;

  typedef   std::map<int,EdgeVector>::iterator        iterator;
  typedef   std::map<int,EdgeVector>::const_iterator  const_iterator;

public:

  HNStructureQ1() {};
  virtual ~HNStructureQ1() {}
  virtual void SparseStructureDiag(SparseStructure* S) const=0;
  virtual void ReInit(const MeshInterface* m)=0;
  
  virtual void MatrixDiag(int ncomp, MatrixInterface& A) const=0;
  virtual void Average(Vector& u) const=0;
  virtual void Distribute(Vector& u) const=0;
  virtual void Zero(Vector& u) const=0;
  virtual bool ZeroCheck(const Vector& u) const=0;

  virtual void CondenseHanging(nvector<int>& indices) const=0;
  virtual void CondenseHanging(EntryMatrix& E, nvector<int>& indices) const=0;
  virtual void CondenseHangingPatch(EntryMatrix& E, nvector<int>& indices) const=0;
};

/*-----------------------------------------*/


#endif
