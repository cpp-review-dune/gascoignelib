#ifndef  __boundarymanager_h
#define  __boundarymanager_h

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

#include  <map>
#include  <string>
#include  "stlio.h"
#include  "gascoigne.h"
#include  "paramfile.h"

/*---------------------------------------------------------------*/

class BoundaryManager
{
 public:

  typedef  std::map<int,Gascoigne::IntVector>::const_iterator  const_iterator;

 protected:

  std::set<int>                 coldir, colneu;
  std::map<int,Gascoigne::IntVector>       dirvec;

  void AddDirichlet(int col, int c)    
    {
      coldir.insert(col);
      dirvec[col].push_back(c);
    }
  void AddNeumann(int col)    
    {
      colneu.insert(col);
    }

 public:

  BoundaryManager() {}
  BoundaryManager(const Gascoigne::ParamFile* pf);

  virtual std::string GetName() const {return "Std";}

  std::ostream& print(std::ostream& s) const;

  virtual const Gascoigne::IntSet&    GetNeumannColors      (     ) const { return colneu;}
  virtual const Gascoigne::IntSet&    GetDirichletColors    (     ) const { return coldir;}
  virtual const Gascoigne::IntVector& GetDirichletComponents(int c) const 
    { 
      std::map<int,Gascoigne::IntVector>::const_iterator p = dirvec.find(c);
      if(p==dirvec.end())
	{
	  std::cerr << "BoundaryManager::Components()\n";
	  std::cerr << "no such color " << c <<std::endl;
	  std::cerr << "dirvec = " << dirvec << std::endl;
	  abort();
	}
      return p->second;
    }
};

#endif
