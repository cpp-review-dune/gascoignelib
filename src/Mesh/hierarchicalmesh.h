#ifndef __hierarchicalmesh_h
#define __hierarchicalmesh_h

#include  <utility>

#include  "triple.h"
#include  <string>
#include  "edge.h" 
#include  "hanglist.h" 
#include  "curvedshapes.h"
#include  "meshinterface.h"
#include  "vertex.h"
#include  "compvector.h"
#include  "gascoigne.h"
#include  "paramfile.h"

/*---------------------------------------------------*/

class HierarchicalMesh : public MeshInterface
{
 protected :

  /*  typedef  */

  typedef  std::pair<int,int>                 pint;
  typedef  triple<int,int,int>                tint;
  typedef  fixarray<2,int>                    EdgeVector;
  typedef  fixarray<4,int>                    FaceVector;
  typedef  std::vector<Edge>                  EdgeVec;
  typedef  Gascoigne::IntSet::iterator        IntSetIt;   
  typedef  Gascoigne::IntSet::const_iterator  IntSetCIt;
 
  /*  Data  */

  int       mnlevels, pdepth, etapatcher;
  Gascoigne::IntVector    vo2n, eo2n, co2n;
  EdgeVec   edges;
  
  void  update_edges(Gascoigne::IntVector&);
  virtual int   FindPatchDepth() const=0;
  virtual void  FillVertexLevels(Gascoigne::IntVector& dst) const=0;
  virtual void  RefineCoarseNodes(Gascoigne::IntSet& dst, const Gascoigne::IntVector& refnodes,
				  const Gascoigne::IntVector& vertexlevel) const=0;
  virtual void  VertexToCells(Gascoigne::IntVector& dst, const Gascoigne::IntSet& src, 
			      const Gascoigne::IntVector& vertexlevel) const=0;
  virtual void VertexToCellsCoarsening(Gascoigne::IntVector& dst, const Gascoigne::IntSet& src, 
			       const Gascoigne::IntVector& vertexlevel) const=0;

  public:

  virtual ~HierarchicalMesh();
  
  int  withfaces;

  HierarchicalMesh();
  HierarchicalMesh(const HierarchicalMesh&);
  HierarchicalMesh& operator=(const HierarchicalMesh&);

  /*  Zugriff  */

  virtual int  nnodes ()        const =0;
  int  nlevels()                const { return 1+mnlevels;}
  int  nedges ()                const { return edges.size();}

  const Gascoigne::IntVector&  Vertexo2n()    const { return vo2n;}
  const Gascoigne::IntVector&  Edgeo2n  ()    const { return eo2n;}
  const Gascoigne::IntVector&  Cello2n  ()    const { return co2n;}

  int Vertexo2n(int i)          const { assert(i<vo2n.size());return vo2n[i];}
  int Edgeo2n  (int i)          const { assert(i<eo2n.size());return eo2n[i];}
  int Cello2n  (int i)          const { assert(i<co2n.size());return co2n[i];}

  const Edge&    edge(int i)    const {assert(i<edges.size());return edges[i];}
  const EdgeVec& edge()         const { return edges;}

  virtual int  level(int i)   const =0;
  virtual bool sleep(int i)   const =0;

  virtual int    Vater(const int i) const
    {assert(0);}
  virtual Gascoigne::IntVector Nachkommen(const int i) const
    {assert(0);}
  virtual Gascoigne::IntVector Geschwister(const int i) const
    {assert (0);} ;
  virtual Gascoigne::IntVector Kinder(const int i) const
    {assert(0);}
    
  void SetParameters(std::string gridname, int patchdepth, int epatcher);
  void ReadFile(const std::string& gridname);
  void ReadParameters(const Gascoigne::ParamFile* pf);
  void global_refine  (int k);
  void random_refine  (double, int k=1);
  void random_patch_refine  (double, int k=1);
  void random_double_patch_refine  (double, int k=1);
  void clear_transfer_lists();
  virtual void   write_gup(const std::string&) const =0;

  virtual int dimension()       const { return 0;}
  virtual int ncells ()         const =0;
  int patchdepth()              const {return pdepth;}
  virtual int nactivedescendants(int i)      const=0;
  virtual Gascoigne::IntVector GetVertices(int c) const=0;
  
  bool CellIsCurved(int iq) const  { return GetBoundaryCellOfCurved(iq)!=-1;}

  virtual std::set<int> GetColors() const=0;

  virtual void read_inp (const std::string&)=0;
  virtual void read_gup (const std::string&)=0;
  virtual void refine(const Gascoigne::IntVector&, const Gascoigne::IntVector&)=0;
  virtual void patch_refine(Gascoigne::IntVector&, Gascoigne::IntVector&)=0;
  virtual void vertex_patch_refine(Gascoigne::IntVector& ref, Gascoigne::IntVector& coarse);
  virtual void vertex_patch_refine(Gascoigne::IntVector&);
  virtual void GetAwakePatchs(std::set<int>&) const =0;
  virtual void GetAwakeCells(std::set<int>&) const =0;
  virtual void ConstructQ2PatchMesh(Gascoigne::IntVector& pm) const=0;
  virtual std::set<int> CellNeighbours(int i) const 
    { std::cerr << "no CellNeighbours"; abort();}

  virtual int  GetBoundaryCellOfCurved(int iq) const { return -1;}

  virtual void AddShape(int col, BoundaryFunction<2>* f) {assert(0);}
  virtual void AddShape(int col, BoundaryFunction<3>* f) {assert(0);}
};

/*---------------------------------------------------*/

#endif
