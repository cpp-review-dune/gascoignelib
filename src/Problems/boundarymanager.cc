#include  "boundarymanager.h"
#include  "filescanner.h"
#include  "stlio.h"

/*-------------------------------------------------------*/

BoundaryManager::BoundaryManager(const ParamFile* pf)
{
  DataFormatHandler DF;
  DF.insert("dirichlet",&coldir);
  DF.insert("neumann"  ,&colneu);
  DF.insert("dirichletcomp",&dirvec);
  FileScanner FS(DF,pf,"BoundaryManager");

  for(const_iterator p=dirvec.begin();p!=dirvec.end();p++)
    {
      if( coldir.find(p->first) == coldir.end() )
	{
	  std::cerr << "BoundaryManager::BoundaryManager()\n";
	  std::cerr << "problem in component data\n";
	  std::cerr << "color not found: " << p->first << std::endl;
	  abort();
	}
    }
  //print(std::cerr);
}

/*-------------------------------------------------------*/

std::ostream& BoundaryManager::print(std::ostream& s) const
{
  s << "BoundaryManager:\n";
  s << "coldir:\t" << coldir << std::endl;
  s << "colneu:\t" << colneu << std::endl;
  s << "dirvec:\t" << dirvec << std::endl;
  return s;
}
