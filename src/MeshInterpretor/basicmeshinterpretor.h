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

class BasicMeshInterpretor : public MeshInterpretorInterface
{
 private:
  
  const MeshInterface*           __MP;
  FemInterface*         __FEM;
  IntegratorInterface*  __INT;
  
 protected:

  mutable EntryMatrix __E;
  mutable Gascoigne::LocalVector __F;
  mutable Gascoigne::LocalVector __U;

  mutable Gascoigne::LocalData   __Q;
  
  const MeshInterface* GetMesh() const {return __MP;}
  const FemInterface* GetFem() const {return __FEM;}
  const IntegratorInterface* GetIntegrator() const {return __INT;}

  IntegratorInterface*& GetIntegratorPointer() {return __INT;}
  FemInterface*& GetFemPointer() {return __FEM;}

  virtual void GlobalToLocal(Gascoigne::LocalVector& U, const Gascoigne::GlobalVector& u, int iq) const {
    GlobalToLocalSingle(U,u,iq);
    GlobalToLocalData(iq);
  }
  virtual void GlobalToLocalData(int iq) const;

  virtual void GlobalToLocalSingle(Gascoigne::LocalVector& U, const Gascoigne::GlobalVector& u, int iq) const;
  virtual void LocalToGlobal(Gascoigne::GlobalVector& f, const Gascoigne::LocalVector& F, int iq, double s) const;
  virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

  virtual nvector<int> GetLocalIndices(int iq) const=0;

 public:
  
  //
  ////  Constructor 
  //
  
  BasicMeshInterpretor();
  ~BasicMeshInterpretor();
  
  void BasicInit(const Gascoigne::ParamFile* pf) {}
  void ReInit   (const MeshInterface* MP) {__MP=MP;}
};


#endif
