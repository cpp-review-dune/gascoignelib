#ifndef __boundaryindexhandler_h
#define __boundaryindexhandler_h

#include <map>
#include <set>
#include  "gascoigne.h"

using namespace Gascoigne;

/*--------------------------------------------------------------*/

class BoundaryIndexHandler
{
 protected:

  typedef std::map<int,IntVector> VecMap;

  IntSet  AllColors;
  VecMap  verteces, cells, localci;

 public:

  void CopySetToVector(const std::vector<IntSet>&,
		       const IntVector&, VecMap&) const;

  void clear();

  const IntSet& GetColors() const {return AllColors;}
  const VecMap& GetVertex() const {return verteces;}
  const VecMap& GetCell()   const {return cells;}
  const VecMap& GetLocal() const  {return localci;}

  IntSet& GetColors()  {return AllColors;}
  VecMap& GetVertex()  {return verteces;}
  VecMap& GetCell()    {return cells;}
  VecMap& GetLocal()   {return localci;}

  const IntVector& Verteces(int col) const;
  const IntVector& Cells   (int col) const;
  const IntVector& Localind(int col) const;

  void Equal(const IntSet& col, const VecMap& v, const VecMap& c, const VecMap& l);
  void check() const;
  friend std::ostream& operator<<(std::ostream &s, const BoundaryIndexHandler& A);
};


/*--------------------------------------------------------------*/

#endif
