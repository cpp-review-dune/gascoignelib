#ifndef __hangcontainer2d_h
#define __hangcontainer2d_h

#include  "hanglist.h"
#include  "gascoigne.h"

using namespace Gascoigne;

/*********************************************************************/

class HangContainer2d
{
 protected:

  typedef  fixarray<2,int>               EdgeVector;
  
  HangList<2>   VertexToBeDeleted;  // ehemals haengende Knoten
  HangList<2>   VertexToBeCreated;  // neue haengende Knoten
  HangList<2>   NotAnyMoreHanging;  // Knoten bleibt aber haengt nicht mehr

  HangList<2>&  Hanging;  // dies sind die aktuellen haengenden Knoten

 public:

  HangContainer2d(HangList<2>& lh) : Hanging(lh) {}

  // Zugriff

  int  NToBeDeleted()   const { return VertexToBeDeleted.size(); }
  int  NToBeCreated()   const { return VertexToBeCreated.size(); }

  const HangList<2>& Deleting()   const { return VertexToBeDeleted;}
  const HangList<2>& Creating()   const { return VertexToBeCreated;}
  const HangList<2>& NotAnyMore() const { return NotAnyMoreHanging;}
        HangList<2>& NotAnyMore()       { return NotAnyMoreHanging;}

  bool ToBeDeleted(const EdgeVector& v) const;
  bool ToBeCreated(const EdgeVector& v) const;

  void make_consistent() { VertexToBeCreated.make_consistent(VertexToBeDeleted);}

  void load_elimination(IntVector&) const;

  void update_olds (const IntVector&, const IntVector&);
  void update_news (const IntVector&,int);

  int  vertex_index (const EdgeVector&) const;

  void ghost_coarse (EdgeVector&,  int, int);
  void ghost_refine (EdgeVector&,  int);

  void NeighbourSwapper();
};

/*********************************************************************/

#endif
