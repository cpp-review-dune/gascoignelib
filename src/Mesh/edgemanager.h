#ifndef __edgemanager_h
#define __edgemanager_h

#include  "edge.h"
#include  "quadlawandorder.h" 
#include  "hangcontainer2d.h" 
#include  "gascoigne.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class EdgeManager
{
  typedef fixarray<2,int>  EdgeVector;

 protected:

  std::vector<Edge>&          edges;
  std::vector<Quad>&          quads;
  const IntVector&            co2n;
        IntVector&            eo2n;

  IntVector                SwappedEdge;
  QuadLawAndOrder          QuadLaO;

  void Update      ();
  void InnerEdges  (const IntSet& CellRefList);
  void OuterEdges  (const HangContainer2d& hangset);
  void OldHangings (HangContainer2d& hangset, const IntSet& CellRefList);
  void SwappedEdges();
  void NeighbourTester() const;

  void BSETest() const;

 public:

  EdgeManager(std::vector<Edge>&, std::vector<Quad>&, const IntVector& con, IntVector& eon);

  const Quad&  quad(int i)           const { return quads[i];}
        Quad&  quad(int i)                 { return quads[i];}

  fixarray<2,int> ChildrenOfEdge(int e) const;

  bool EdgeIsHanging(int e) const;
  bool EdgeIsHanging(const Edge& e) const;

  void LoadEdgeElimination(IntVector& edel, 
			   const IntSet& CellCoarseList,
			   const HangContainer2d& hangset) const;
  void Build( const IntSet& CellRefList, HangContainer2d&);
  void DeleteEdges();
  void InitEdges();
  void SortHangings();
};
}

/*---------------------------------------------------*/

#endif
