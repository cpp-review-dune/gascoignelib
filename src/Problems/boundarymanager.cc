#include  "boundarymanager.h"
#include  "filescanner.h"
#include  "stlio.h"


using namespace std;

/*-------------------------------------------------------*/

namespace Gascoigne
{
BoundaryManager::BoundaryManager(const ParamFile* pf)
{
  cerr << "------------------- BoundaryManager::BoundaryManager(const ParamFile* pf) depreciated!\n";
  BasicInit(pf);
}

/*-------------------------------------------------------*/

void BoundaryManager::BasicInit(const ParamFile* pf)
{
  DataFormatHandler DF;
  DF.insert("dirichlet",&coldir);
  DF.insert("neumann"  ,&colneu);
  DF.insert("robin"    ,&colrob);
  DF.insert("dirichletcomp",&dirvec);
  FileScanner FS(DF,pf,"BoundaryManager");

  for(const_iterator p=dirvec.begin();p!=dirvec.end();p++)
    {
      if( coldir.find(p->first) == coldir.end() )
	{
	  cerr << "BoundaryManager::BoundaryManager()\n";
	  cerr << "problem in component data\n";
	  cerr << "color not found: " << p->first << endl;
	  abort();
	}
    }
  //print(cerr);
}

/*-------------------------------------------------------*/

ostream& BoundaryManager::print(ostream& s) const
{
  s << "BoundaryManager:\n";
  s << "coldir:\t" << coldir << endl;
  s << "colneu:\t" << colneu << endl;
  s << "dirvec:\t" << dirvec << endl;
  return s;
}
}
