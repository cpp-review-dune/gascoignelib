#ifndef  __MeshAgent_h
#define  __MeshAgent_h

#include  "meshagentinterface.h"
#include  "hierarchicalmesh.h"
#include  "gascoignemultigridmesh.h"
#include  "boundaryfunction.h"

/*-----------------------------------------*/

class MeshAgent : public virtual MeshAgentInterface
{
protected:

  int _dimension;

  HierarchicalMesh*         HMP;
  GascoigneMultiGridMesh*   GMG;

  virtual GascoigneMultiGridMesh* NewMultiGridMesh() {return new GascoigneMultiGridMesh;}

  virtual void ReInit();

  GascoigneMesh*  GMesh(int l) { return GMG->GetGascoigneMesh(l);}

public:
    
  MeshAgent();
  ~MeshAgent();

  void BasicInit(const Gascoigne::ParamFile* pf);
  void BasicInit(int dim, std::string inpname, int prerefine);
  void BasicInit(std::string inpname, int prerefine, std::map<int,BoundaryFunction<2>* >& curved);
  void BasicInit(std::string inpname, int prerefine, std::map<int,BoundaryFunction<3>* >& curved);

  const GascoigneMultiGridMesh& GetMultiGrid() const {return *GMG;}
  GascoigneMultiGridMesh& GetMultiGrid() {return *GMG;}

  int nnodes() const {return GMG->GetGascoigneMesh(0)->nnodes();}
  int ncells() const {return GMG->GetGascoigneMesh(0)->ncells();}
  int nlevels() const {return GMG->nlevels();}

  const MeshInterface* GetMesh()      const { return GMG->GetGascoigneMesh(0);}
  const MeshInterface* GetMesh(int l) const { return GMG->GetGascoigneMesh(l);}

  void read_gup(const std::string& fname);
  void write_gup(const std::string& fname) const;
  void global_refine(int n);
  void random_patch_refine(double p, int n);
  void refine_nodes(nvector<int>& refnodes, nvector<int>& coarsenodes);
  void refine_nodes(nvector<int>& refnodes);
  void refine_cells(nvector<int>& ref);

  const GascoigneMeshTransfer* GetTransfer(int l) const {return GMG->GetTransfer(l);}

  const HierarchicalMesh* getHierarchicalMesh() const { return HMP;}
};


#endif
