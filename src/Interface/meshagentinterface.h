#ifndef  __MeshAgentInterface_h
#define  __MeshAgentInterface_h


/////////////////////////////////////////////
////
////@brief
////  ... comments MeshAgentInterface

////
////
/////////////////////////////////////////////

#include  "meshinterface.h"
#include  "multigridmeshinterface.h"
#include  "meshtransferinterface.h"
#include  "paramfile.h"
#include  <string>

namespace Gascoigne
{
  class MeshAgentInterface
  {
    private:

    protected:

    public:
      //
      ////  Con(De)structor 
      //
      MeshAgentInterface() {}
      virtual ~MeshAgentInterface() {}

      virtual void BasicInit(const ParamFile* pf)=0;

      virtual void SetDefaultValues(int dim, std::string gridname, int prerefine)=0;

      virtual int nnodes() const=0;
      virtual int ncells() const=0;
      virtual int nlevels() const=0;

      virtual void write_gup(const std::string& fname) const=0;
      virtual const MeshInterface* GetMesh(int l) const=0;

      virtual void global_patch_coarsen(int n)=0;
      virtual void global_refine(int n)=0;
      virtual void refine_nodes(IntVector& refnodes, IntVector& coarsenodes)=0;
      virtual void random_patch_refine(double p, int n)=0;
      virtual void random_patch_coarsen(double p, int n)=0;
      virtual const MeshTransferInterface* GetTransfer(int l) const=0;
  };
}

#endif
