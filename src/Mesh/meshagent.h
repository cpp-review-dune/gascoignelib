#ifndef  __MeshAgent_h
#define  __MeshAgent_h

#include  "meshagentinterface.h"
#include  "hierarchicalmesh.h"
#include  "gascoignemultigridmesh.h"
#include  "boundaryfunction.h"
#include  "stdperiodicmapping.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
    typedef std::map<int,int> IntMap;

class MeshAgent : public virtual MeshAgentInterface
{
private:

  map<int,BoundaryFunction<2>* > _curved2d;
  map<int,BoundaryFunction<3>* > _curved3d;

  map<int,map<int,PeriodicMapping*> > _periodicMaps;
  IntVector _periodicCols;

  //Fuer die Zuordnung GM Nr auf altem Gitter zu GM Nr. auf neuem Gitter
  IntVector _cl2g, _celll2g;
  IntVector _fathers;//GM Nr zu HM nr.
  IntMap _cg2l, _cellg2l;
  map<int,set<int> > _co2n;
  bool _goc2nc;
  nvector<IntVector> _q4patch,_q4toq2;

protected:

  int GetDimension() const {return HMP->dimension();}	
	
  HierarchicalMesh*         HMP;
  GascoigneMultiGridMesh*   GMG;

  virtual GascoigneMultiGridMesh* NewMultiGridMesh() {return new GascoigneMultiGridMesh;}

  virtual void ReInit();
  virtual void BuildQ4PatchList(const IntVector &patchl2g);
  virtual void AssemblePeriodicBoundaries();

  GascoigneMesh*  GMesh(int l) { return GMG->GetGascoigneMesh(l);}

public:
    
  MeshAgent();
  ~MeshAgent();

  void AddShape(int col, BoundaryFunction<2>* f) { _curved2d[col] = f;}
  void AddShape(int col, BoundaryFunction<3>* f) { _curved3d[col] = f;}

        map<int,BoundaryFunction<2>* >& GetShapes2d()       { return _curved2d; }
        map<int,BoundaryFunction<3>* >& GetShapes3d()       { return _curved3d; }
  const map<int,BoundaryFunction<2>* >& GetShapes2d() const { return _curved2d; }
  const map<int,BoundaryFunction<3>* >& GetShapes3d() const { return _curved3d; }

  void AddPeriodicMapping(int col, int col2, PeriodicMapping* map) { _periodicMaps[col][col2] = map; }

  void BasicInit(const ParamFile* pf);
  void BasicInit(const std::string& gridname, int dim, int patchdepth, int epatcher, bool goc2nc = false);

  const GascoigneMultiGridMesh& GetMultiGrid() const {return *GMG;}
  GascoigneMultiGridMesh& GetMultiGrid() {return *GMG;}

        HierarchicalMesh* GetHierarchicalMesh()       { return HMP;}
  const HierarchicalMesh* GetHierarchicalMesh() const { return HMP;}

  int nnodes() const {return GMG->GetGascoigneMesh(0)->nnodes();}
  int ncells() const {return GMG->GetGascoigneMesh(0)->ncells();}
  int nlevels() const {return GMG->nlevels();}

  const MeshInterface* GetMesh()      const { return GMG->GetGascoigneMesh(0);}
  const MeshInterface* GetMesh(int l) const { return GMG->GetGascoigneMesh(l);}

  void read_gup(const std::string& fname);
  void read_gip(const std::string& fname);
  void write_gup(const std::string& fname) const;
  void write_gip(const std::string& fname) const;
  void write_inp(const std::string& fname) const;
  void global_refine(int n);
  void global_patch_coarsen(int n);
  void random_patch_coarsen(double p, int n);
  void random_patch_refine(double p, int n);
  void refine_nodes(IntVector& refnodes, IntVector& coarsenodes);
  void refine_nodes(IntVector& refnodes);
  void refine_cells(IntVector& ref);

  const GascoigneMeshTransfer* GetTransfer(int l) const {return GMG->GetTransfer(l);}

  const std::set<int> Cello2n(int i)const;
  const int Cello2nFather(int i)const;
  void ClearCl2g() { _cl2g.clear(); }
  const bool Goc2nc() const{ return _goc2nc;}

  const IntVector& Celll2g() const { return _celll2g; }
  const IntMap& Cellg2l()    const { return _cellg2l; }
};
}

#endif
