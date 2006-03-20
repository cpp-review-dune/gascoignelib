#ifndef  __MeshAgentInterface_h
#define  __MeshAgentInterface_h


#include  "meshinterface.h"
#include  "multigridmeshinterface.h"
#include  "meshtransferinterface.h"
#include  "paramfile.h"
#include  <string>

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments MeshAgentInterface

  ////
  ////
  /////////////////////////////////////////////

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

      virtual int GetDimension() const=0;

      virtual void BasicInit(const ParamFile* pf)=0;
      virtual void BasicInit(const std::string& gridname, int dim, int patchdepth, int epatcher, bool goc2nc=false)=0;

      virtual int nnodes() const=0;
      virtual int ncells() const=0;
      virtual int nlevels() const=0;

      virtual void read_gup(const std::string& fname)=0;
      virtual void read_gip(const std::string& fname)=0;
      virtual void write_gup(const std::string& fname) const=0;
      virtual void write_gip(const std::string& fname) const=0;
      virtual void write_inp(const std::string& fname) const=0;
      virtual const MeshInterface* GetMesh(int l) const=0;

      virtual void global_patch_coarsen(int n)=0;
      virtual void global_refine(int n)=0;
      virtual void refine_nodes(IntVector& refnodes, IntVector& coarsenodes)=0;
      virtual void refine_nodes(IntVector& refnodes)=0;
      virtual void refine_cells(IntVector& refnodes)=0;
      virtual void random_patch_refine(double p, int n)=0;
      virtual void random_patch_coarsen(double p, int n)=0;
      virtual const MeshTransferInterface* GetTransfer(int l) const=0; 
      virtual const std::set<int> Cello2n(int i)const=0;
      virtual const int Cello2nFather(int i)const=0;
      virtual const bool Goc2nc()const=0;
  };
}

#endif
