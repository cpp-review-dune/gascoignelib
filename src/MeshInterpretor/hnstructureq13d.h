#ifndef  __HNStructureQ13d_h
#define  __HNStructureQ13d_h

#include  "hnstructureq12d.h"
#include  "entrymatrix.h"

/*-----------------------------------------*/

class HNStructureQ13d : public HNStructureQ12d
{
protected:

  typedef   fixarray<9,int>    FaceVector;

  typedef   std::map<int,FaceVector>::iterator        fiterator;
  typedef   std::map<int,FaceVector>::const_iterator  const_fiterator;

  const std::map<int,FaceVector>*        faces;

  void CondenseHanging2er(nvector<int>& indices) const;
  void CondenseHanging4er(nvector<int>& indices) const;

  void CondenseHanging2er(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHanging4er(EntryMatrix& E, nvector<int>& indices) const;

  fixarray<4,int> GetHangingFace(int i) const;
  fixarray<2,int> GetHangingEdge(int i) const;

public:

  ~HNStructureQ13d();
  HNStructureQ13d();

  int nhnodes() const { return edges->size() + faces->size();} 

  void ReInit(const MeshInterface* m);
  int   hanging(int i) const;

  void MatrixDiag(int ncomp, MatrixInterface& A) const;
  void SparseStructureDiag(SparseStructure* A) const;
  
  void Average(Gascoigne::GlobalVector& u) const;
  void Distribute(Gascoigne::GlobalVector& u) const;
  void Zero(Gascoigne::GlobalVector& u) const;
  bool ZeroCheck(const Gascoigne::GlobalVector& u) const;
  
  void Couplings(nvector<int>& indices) const;
  
  void CondenseHanging(nvector<int>& indices) const;
  void CondenseHanging(EntryMatrix&, nvector<int>&) const;
  void CondenseHangingPatch(EntryMatrix& E, nvector<int>& indices) const {}
};

#endif
