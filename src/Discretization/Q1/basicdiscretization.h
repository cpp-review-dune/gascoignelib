#ifndef  __BasicDiscretization_h
#define  __BasicDiscretization_h


#include  "discretizationinterface.h"
#include  "feminterface.h"
#include  "integratorinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments BasicDiscretization

///
///
/////////////////////////////////////////////

class BasicDiscretization : public DiscretizationInterface
{
 private:
   const MeshInterface*  __MP;
   mutable GlobalData    __q;
  
 protected:
   mutable EntryMatrix __E;
   mutable LocalVector __F;
   mutable LocalVector __U;

   mutable LocalNodeData        __QN;
   mutable LocalParameterData   __QP;
   mutable LocalCellData        __QC;
   
   virtual const GlobalData& GetGlobalData() const {return __q;}
   virtual void SetGlobalData(const GlobalData& q) const {__q = q;}

   virtual const MeshInterface* GetMesh() const { assert(__MP); return __MP;}

   virtual void GlobalToGlobalData() const;
   virtual void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const {
     GlobalToLocalSingle(U,u,iq);
     GlobalToLocalData(iq);
   }
   virtual void GlobalToLocalData(int iq) const;
   virtual void GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const;
   virtual void GlobalToLocalCell(LocalCellVector& U, const GlobalCellVector& u, int iq) const;

   virtual void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const;
   virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

   virtual IntVector GetLocalIndices(int iq) const=0;

 public:
   //
   ////  Constructor 
   //
   BasicDiscretization();
   ~BasicDiscretization();
   
   void BasicInit(const ParamFile* pf) {}
   void ReInit   (const MeshInterface* MP) {__MP=MP;}
   
   virtual void AddNodeVector(const std::string& name, const GlobalVector* q) const {
     __q.AddNodeVector(name,q);
   }
   virtual void DeleteNodeVector(const std::string& name) const {
     __q.DeleteNodeVector(name);
   }
   virtual void AddCellVector(const std::string& name, const GlobalCellVector* q) const {
     __q.AddCellVector(name,q);
   }
   virtual void DeleteCellVector(const std::string& name) const {
     __q.DeleteCellVector(name);
   }
   virtual void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const {
     __q.AddParameterVector(name,q);
   }
   virtual void DeleteParameterVector(const std::string& name) const {
     __q.DeleteParameterVector(name);
   }

   void HNAverageData() const;
   void HNZeroData   () const;
};
}

#endif
