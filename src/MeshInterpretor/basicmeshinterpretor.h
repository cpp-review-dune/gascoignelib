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

using namespace Gascoigne;

class BasicMeshInterpretor : public MeshInterpretorInterface
{
 private:
  
  const MeshInterface*           __MP;
  FemInterface*         __FEM;
  IntegratorInterface*  __INT;
  
 protected:

  mutable EntryMatrix __E;
  mutable LocalVector __F;
  mutable LocalVector __U;

  mutable LocalData   __Q;
  
  const MeshInterface* GetMesh() const {return __MP;}
  const FemInterface* GetFem() const {return __FEM;}
  const IntegratorInterface* GetIntegrator() const {return __INT;}

  IntegratorInterface*& GetIntegratorPointer() {return __INT;}
  FemInterface*& GetFemPointer() {return __FEM;}

  virtual void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const {
    GlobalToLocalSingle(U,u,iq);
    GlobalToLocalData(iq);
  }
  virtual void GlobalToLocalData(int iq) const;

  virtual void GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const;
  virtual void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const;
  virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

  virtual nvector<int> GetLocalIndices(int iq) const=0;

 public:
  
  //
  ////  Constructor 
  //
  
  BasicMeshInterpretor();
  ~BasicMeshInterpretor();
  
  void BasicInit(const std::string& paramfile);
  void ReInit   (const MeshInterface* MP) {__MP=MP;}
};


#endif
