#ifndef __edgemanager_h
#define __edgemanager_h

#include  <set>
#include  "edge.h"
#include  "quadlawandorder.h" 
#include  "hangcontainer2d.h" 
#include  "gascoigne.h"

/*---------------------------------------------------*/

class EdgeManager
{
  typedef fixarray<2,int>  EdgeVector;

 protected:

  std::vector<Edge>&          edges;
  std::vector<Quad>&          quads;
  const Gascoigne::IntVector&            co2n;
        Gascoigne::IntVector&            eo2n;

  Gascoigne::IntVector                SwappedEdge;
  QuadLawAndOrder          QuadLaO;

  void Update      ();
  void InnerEdges  (const Gascoigne::IntSet& CellRefList);
  void OuterEdges  (const HangContainer2d& hangset);
  void OldHangings (HangContainer2d& hangset, const Gascoigne::IntSet& CellRefList);
  void SwappedEdges();
  void NeighbourTester() const;

  void BSETest() const;

 public:

  EdgeManager(std::vector<Edge>&, std::vector<Quad>&, const Gascoigne::IntVector& con, Gascoigne::IntVector& eon);

  const Quad&  quad(int i)           const { return quads[i];}
        Quad&  quad(int i)                 { return quads[i];}

  fixarray<2,int> ChildrenOfEdge(int e) const;

  bool EdgeIsHanging(int e) const;
  bool EdgeIsHanging(const Edge& e) const;

  void LoadEdgeElimination(Gascoigne::IntVector& edel, 
			   const Gascoigne::IntSet& CellCoarseList,
			   const HangContainer2d& hangset) const;
  void Build( const Gascoigne::IntSet& CellRefList, HangContainer2d&);
  void DeleteEdges();
  void InitEdges();
  void SortHangings();
};

/*---------------------------------------------------*/

#endif
