#ifndef  __BasicMeshInterpretor_h
#define  __BasicMeshInterpretor_h

/////////////////////////////////////////////
///
///@brief
///  ... comments BasicMeshInterpretor

///
///
/////////////////////////////////////////////


#include  "meshinterpretorinterface.h"
#include  "feminterface.h"
#include  "integratorinterface.h"

namespace Gascoigne
{
class BasicMeshInterpretor : public MeshInterpretorInterface
{
 private:
  
  const MeshInterface*  __MP;
  
 protected:

  mutable EntryMatrix            __E;
  mutable LocalVector __F;
  mutable LocalVector __U;

  mutable LocalNodeData        __Q;
  mutable LocalParameterData   __qq;
  
  const MeshInterface* GetMesh() const { assert(__MP); return __MP;}

  virtual void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const {
    GlobalToLocalSingle(U,u,iq);
    GlobalToLocalData(iq);
  }
  virtual void GlobalToLocalData(int iq) const;
  virtual void GlobalToGlobalData() const;

  virtual void GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const;
  virtual void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const;
  virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

  virtual IntVector GetLocalIndices(int iq) const=0;

 public:
  
  //
  ////  Constructor 
  //
  
  BasicMeshInterpretor();
  ~BasicMeshInterpretor();
  
  void BasicInit(const ParamFile* pf) {}
  void ReInit   (const MeshInterface* MP) {__MP=MP;}
  
  void HNAverageData() const;
  void HNZeroData   () const;
};
}

#endif
