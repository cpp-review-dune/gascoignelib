#ifndef  __MeshAgent_h
#define  __MeshAgent_h

#include  "meshagentinterface.h"
#include  "hierarchicalmesh.h"
#include  "gascoignemultigridmesh.h"
#include  "boundaryfunction.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
class MeshAgent : public virtual MeshAgentInterface
{
private:

  int _dimension,_prerefine;
  std::string    _gridname;
  map<int,BoundaryFunction<2>* > _curved2d;
  map<int,BoundaryFunction<3>* > _curved3d;

protected:

		int GetDimension() const {return _dimension;}	
	int GetPrerefine() const {return _prerefine;}	
	std::string GetGridname() const {return _gridname;}

		
  HierarchicalMesh*         HMP;
  GascoigneMultiGridMesh*   GMG;

  virtual GascoigneMultiGridMesh* NewMultiGridMesh() {return new GascoigneMultiGridMesh;}

  virtual void ReInit();
  virtual void ReadMesh(int dim, std::string meshname, int prerefine);

  GascoigneMesh*  GMesh(int l) { return GMG->GetGascoigneMesh(l);}

public:
    
  MeshAgent();
  ~MeshAgent();

  void AddShape(int col, BoundaryFunction<2>* f) { _curved2d[col] = f;}
  void AddShape(int col, BoundaryFunction<3>* f) { _curved3d[col] = f;}

  void ReadParamFile(const ParamFile* paramfile);
  void BasicInit(const ParamFile* pf);

  void SetDefaultValues(int dimension, std::string gridname, int prerefine);

  const GascoigneMultiGridMesh& GetMultiGrid() const {return *GMG;}
  GascoigneMultiGridMesh& GetMultiGrid() {return *GMG;}

  const HierarchicalMesh* GetHierarchicalMesh() const { return HMP;}

  int nnodes() const {return GMG->GetGascoigneMesh(0)->nnodes();}
  int ncells() const {return GMG->GetGascoigneMesh(0)->ncells();}
  int nlevels() const {return GMG->nlevels();}

  const MeshInterface* GetMesh()      const { return GMG->GetGascoigneMesh(0);}
  const MeshInterface* GetMesh(int l) const { return GMG->GetGascoigneMesh(l);}

  void read_gup(const std::string& fname);
  void write_gup(const std::string& fname) const;
  void write_inp(const std::string& fname) const;
  void global_refine(int n);
  void global_patch_coarsen(int n);
  void random_patch_coarsen(double p, int n);
  void random_patch_refine(double p, int n);
  void refine_nodes(IntVector& refnodes, IntVector& coarsenodes);
  void refine_nodes(IntVector& refnodes);
  void refine_cells(IntVector& ref);

  const GascoigneMeshTransfer* GetTransfer(int l) const {return GMG->GetTransfer(l);}
};
}

#endif
