#ifndef  __boundarymanager_h
#define  __boundarymanager_h

#include  <map>
#include  <string>
#include  "stlio.h"
#include  "gascoigne.h"
#include  "paramfile.h"


namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
/// Management of boundary colors

/// According to the input mesh, parts of the boundary are associated 
/// with a number (color). This number is used in different application 
/// classes to identify a certain part of the boundary.
/// The class BoundaryManager administrates a list of numbers for 
/// Dirichlet and Neumann colors.
//
//////////////////////////////////////////////

class BoundaryManager
{
 protected:

  IntSet                   _colsDirichlet, _colsRightHandSide, _colsEquation;
  std::map<int,IntVector>  _compsDirichlet;

 public:

  BoundaryManager() {}
  virtual ~BoundaryManager() {}

  virtual void BasicInit(const ParamFile* pf);

  virtual std::string GetName() const {return "Std";}

  void AddDirichletData(int col, int c)    
    {
      _colsDirichlet.insert(col);
      _compsDirichlet[col].push_back(c);
    }
  void AddBoundaryRightHandSide(int col)    
    {
      _colsRightHandSide.insert(col);
    }
  void AddBoundaryEquation(int col)    
    {
      _colsEquation.insert(col);
    }

  std::ostream& print(std::ostream& s) const;

  virtual const IntSet& GetBoundaryRightHandSideColors() const { return _colsRightHandSide; }
  virtual const IntSet& GetBoundaryEquationColors     () const { return _colsEquation; }
  virtual const IntSet& GetDirichletDataColors        () const { return _colsDirichlet; }

  virtual const IntVector& GetDirichletDataComponents(int c) const 
    { 
      std::map<int,IntVector>::const_iterator p = _compsDirichlet.find(c);
      if(p==_compsDirichlet.end())
	{
	  std::cerr << "BoundaryManager::GetDirichletComponents()" << std::endl;
	  std::cerr << "No such color " << c <<std::endl;
	  std::cerr << "components = " << _compsDirichlet << std::endl;
	  abort();
	}
      return p->second;
    }
};
}

#endif
