#ifndef  __levelmesh2d_h
#define  __levelmesh2d_h

#include  "hierarchicalmesh2d.h"
#include  "index.h"
#include  "boundaryindexhandler.h"

/*--------------------------------------------------------------*/

class LevelMesh2d : public Index
{
 protected:

  typedef std::map<int,fixarray<3,int> >    QuadraticHNStructure3;

  const HierarchicalMesh2d*  HMP;

  void check_leveljump() const;
  void fill_opis(Gascoigne::IntSet& dst, Gascoigne::IntSet& oldquads) const;
  void fill_enkel(Gascoigne::IntSet& dst, const Quad& Q) const;
  void fill_childs(Gascoigne::IntSet& dst, const Quad& Q) const;
  bool EnkelUniform(const Quad& Q) const;
  bool BuildFathers(std::set<int>&  Vaeter) const;
  void InitCells(int n);
  void InitNodes(int n);
  void InitEdges(int n);

 public:
   
   /*----- Constructor -----*/

  LevelMesh2d(const HierarchicalMesh* hmp);
  ~LevelMesh2d();

  const HierarchicalMesh2d* GetHierarchicalMesh() const {return  HMP;}
  
  int             ncells  ()     const  { return Index::QuadSize(); }
  const Quad&     quad    (int i) const { return HMP->quad(Quadl2g(i));}
  const Vertex2d& vertex2d(int i) const { return HMP->vertex2d(Vertexl2g(i));}

  int vertex_of_cell(int i, int j) const 
    { return Vertexg2l(HMP->vertex_of_cell(Quadl2g(i),j)); }

  /*----- Functions -----*/

  bool EdgeIsHangingGlobalIndex(int i) const;

  void BasicInit(const Gascoigne::IntSet& n, const Gascoigne::IntSet& o);

  /*----- Functions for patches -----*/

  void construct_lists(Gascoigne::IntSet& newquads, Gascoigne::IntSet& oldquads) const;
  void ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2) const;
  void InitBoundaryHandler(BoundaryIndexHandler& BI) const;
  void ConstructIndOfPatch(nvector<Gascoigne::IntVector>& dstv) const;
  bool ConstructCellIndOfPatch(nvector<int>& dstc) const;
};

/*---------------------------------------------------*/

#endif
