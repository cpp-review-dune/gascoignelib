#include  "gascoignemultigridmesh.h"
#include  "gascoignemesh2d.h"
#include  "gascoignemesh3d.h"
#include  "gascoignemeshtransfer2d.h"
#include  "gascoignemeshtransfer3d.h"

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMultiGridMesh::GascoigneMultiGridMesh()
{
}

/*-----------------------------------------*/

GascoigneMultiGridMesh::~GascoigneMultiGridMesh()
{
  for(int l=0;l<M.size();l++)  
    {
      if (M[l]!=NULL)  {delete M[l]; M[l]=NULL;}
    }
  for(int l=0;l<T.size();l++)  
    {
      if (T[l]!=NULL)  {delete T[l]; T[l]=NULL;}
    }
}

/*-----------------------------------------*/

GascoigneMesh* GascoigneMultiGridMesh::NewMesh(int dim)
{
  assert(2<=dim && dim<=3);

  if(dim==2)
    return new GascoigneMesh2d;
  else if(dim==3)
    return new GascoigneMesh3d;

  abort();
}

/*-----------------------------------------*/

GascoigneMeshTransfer* GascoigneMultiGridMesh::NewTransfer(int dim)
{
  assert(2<=dim && dim<=3);

  if(dim==2)
    {
      return new GascoigneMeshTransfer2d;
    }
  else if(dim==3)
    {
      return new GascoigneMeshTransfer3d;
    }

  abort();
}

/*-----------------------------------------*/

void GascoigneMultiGridMesh::ReInit(int dim, int nlevels)
{
  // Mesh
  for(int l=0;l<M.size();l++)
    {
      if (M[l]!=NULL) { delete M[l]; M[l]=NULL;}
    }
  M.resize(nlevels);
  for(int l=0;l<M.size();l++)
    {
      M[l] = NewMesh(dim);
    }
  // Transfer
  for(int l=0;l<T.size();l++)
    {
      if(T[l]!=NULL) {delete T[l]; T[l]=NULL;}
    }
  T.resize(nlevels-1);
  for(int l=0;l<T.size();l++)
    {
      T[l] = NewTransfer(dim);
    }
}
}
