#ifndef  __levelmesh3d_h
#define  __levelmesh3d_h

#include  "hierarchicalmesh3d.h"
#include  "index.h"
#include  "boundaryindexhandler.h"
#include  "gascoigne.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{
class LevelMesh3d : public Index
{
 protected:

  typedef std::map<int,fixarray<3,int> >  QuadraticHNStructure3;
  typedef std::map<int,fixarray<9,int> >  QuadraticHNStructure9;

  const HierarchicalMesh3d* HMP;

  void check_leveljump() const;
  void fill_opis(IntSet& dst, IntSet& oldquads) const;
  void fill_enkel(IntSet& dst, const Hex& Q) const;
  void fill_childs(IntSet& dst, const Hex& Q) const;
  bool EnkelUniform(const Hex& Q) const;
  bool BuildFathers(std::set<int>&  Vaeter) const;
  int  hexindex  (int i) const { return Index::Hexl2g(i);}
  void InitCells(int n);
  void InitNodes(int n);
  void InitEdges(int n);

 public:
   
  LevelMesh3d(const HierarchicalMesh* hmp);
  ~LevelMesh3d();

  int ncells   ()         const  { return Index::HexSize(); }
  int dimension()         const  { return HMP->dimension();}

  int vertex_of_cell(int i, int ii) const 
    { return Index::Vertexg2l(HMP->vertex_of_cell(hexindex(i),ii)); }

  const Vertex3d& vertex3d(int i) const 
    { return HMP->vertex3d(Index::Vertexl2g(i)); }

  /*----- Functions -----*/

  const Hex&   hex  (int i) const { return HMP->hex(hexindex(i));}

  void BasicInit(const IntSet& n, const IntSet& o);

  /*----- Functions fuer Patch-----*/

  void construct_lists(IntSet& newquads, IntSet& oldquads) const;

  void ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2,
					  QuadraticHNStructure9& hnq2face) const;

  void InitBoundaryHandler(BoundaryIndexHandler& BI) const;
  bool ConstructCellIndOfPatch(IntVector& dstc) const;
  void ConstructIndOfPatch(nvector<IntVector>& dstv) const;
};
}

/*---------------------------------------------------*/

#endif
