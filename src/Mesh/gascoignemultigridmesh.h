#ifndef  __GascoigneMultiGridMesh_h
#define  __GascoigneMultiGridMesh_h

#include  <vector>
#include  "gascoignemesh.h"
#include  "gascoignemeshtransfer.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMultiGridMesh
{
protected:

  std::vector<GascoigneMesh*>          M;
  std::vector<GascoigneMeshTransfer*>  T;

  virtual GascoigneMesh* NewMesh(int dim);
  virtual GascoigneMeshTransfer* NewTransfer(int dim);

public:
  
  virtual void ReInit(int dim, int nlevels);

  GascoigneMultiGridMesh();
  virtual ~GascoigneMultiGridMesh();

  int nlevels() const {return M.size();}

  const MeshInterface& operator()(int l) const {
    assert((l>=0)&&(l<M.size()));
    return *M[l];}

  const GascoigneMesh* GetGascoigneMesh(int l) const {
    assert((l>=0)&&(l<M.size()));
    return M[l];}
  GascoigneMesh*       GetGascoigneMesh(int l) {
    assert((l>=0)&&(l<M.size()));
    return M[l];}

  GascoigneMeshTransfer* GetTransfer(int l) {
    assert((l>=0)&&(l<T.size()));
    return T[l];}
  const GascoigneMeshTransfer* GetTransfer(int l) const {
    assert((l>=0)&&(l<T.size()));
    return T[l];}
};
}

#endif
