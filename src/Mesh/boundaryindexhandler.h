#ifndef __boundaryindexhandler_h
#define __boundaryindexhandler_h

#include <map>
#include <set>
#include  "gascoigne.h"

/*--------------------------------------------------------------*/

class BoundaryIndexHandler
{
 protected:

  typedef std::map<int,Gascoigne::IntVector> VecMap;

  Gascoigne::IntSet  AllColors;
  VecMap  verteces, cells, localci;

 public:

  void CopySetToVector(const std::vector<Gascoigne::IntSet>&,
		       const Gascoigne::IntVector&, VecMap&) const;

  void clear();

  const Gascoigne::IntSet& GetColors() const {return AllColors;}
  const VecMap& GetVertex() const {return verteces;}
  const VecMap& GetCell()   const {return cells;}
  const VecMap& GetLocal() const  {return localci;}

  Gascoigne::IntSet& GetColors()  {return AllColors;}
  VecMap& GetVertex()  {return verteces;}
  VecMap& GetCell()    {return cells;}
  VecMap& GetLocal()   {return localci;}

  const Gascoigne::IntVector& Verteces(int col) const;
  const Gascoigne::IntVector& Cells   (int col) const;
  const Gascoigne::IntVector& Localind(int col) const;

  void Equal(const Gascoigne::IntSet& col, const VecMap& v, const VecMap& c, const VecMap& l);
  void check() const;
  friend std::ostream& operator<<(std::ostream &s, const BoundaryIndexHandler& A);
};


/*--------------------------------------------------------------*/

#endif
