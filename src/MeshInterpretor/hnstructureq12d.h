#ifndef  __HNStructureQ12d_h
#define  __HNStructureQ12d_h

#include  "hnstructureq1.h"

/*-----------------------------------------*/

class HNStructureQ12d : public virtual HNStructureQ1
{
protected:

  const std::map<int,EdgeVector>*        edges;
  nvector<double>                        wei;
  fixarray<4,fixarray<3,int> >           lnoe, lnop;

  double weight(int i) const { return wei[i];}
  int hanging(int i) const;
  int nhnodes() const {return edges->size();} 
  const EdgeVector& regular_nodes(int i) const;
  
public:

  HNStructureQ12d();
  ~HNStructureQ12d() {}
  void SparseStructureDiag(SparseStructure* S) const;
  void ReInit(const MeshInterface* m);
  
  void MatrixDiag(int ncomp, MatrixInterface& A) const;
  void Average(Vector& u) const;
  void Distribute(Vector& u) const;
  void Zero(Vector& u) const;
  bool ZeroCheck(const Vector& u) const;
  
  void CondenseHanging(nvector<int>& indices) const;
  void CondenseHanging(EntryMatrix& E, nvector<int>& indices) const;
  void CondenseHangingPatch(EntryMatrix& E, nvector<int>& indices) const;
};

#endif
