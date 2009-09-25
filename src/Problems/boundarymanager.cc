#include  "boundarymanager.h"
#include  "filescanner.h"
#include  "stlio.h"


using namespace std;

/*-------------------------------------------------------*/

namespace Gascoigne
{
void BoundaryManager::BasicInit(const ParamFile* pf)
{
  IntSet colsNeumann, colsRobin;
  
  DataFormatHandler DF;
  DF.insert("neumann"      ,&colsNeumann); //not longer supported
  DF.insert("robin"        ,&colsRobin);   //not longer supported

  DF.insert("equation"     ,&_colsEquation);
  DF.insert("righthandside",&_colsRightHandSide);
  DF.insert("functional"   ,&_colsFunctional);
  DF.insert("dirichlet"    ,&_colsDirichlet);
  DF.insert("dirichletcomp",&_compsDirichlet);
  DF.insert("periodic"          ,&_colsPeriodic);
  DF.insert("periodiccomp"      ,&_compsPeriodic);
  FileScanner FS(DF,pf,"BoundaryManager");

  if(colsNeumann.size() || colsRobin.size())
  {
    cerr << "\"neumann\" and \"robin\" are not longer supported!" << endl;
    cerr << "Use \"righthandside\" and \"equation\" instead." << endl;
    abort();
  }

  map<int,IntVector>::const_iterator p = _compsDirichlet.begin();
  for(;p!=_compsDirichlet.end();p++)
    {
      if( _colsDirichlet.find(p->first) == _colsDirichlet.end() )
	{
	  cerr << "BoundaryManager::BoundaryManager()\n";
	  cerr << "problem in Dirichlet component data\n";
	  cerr << "color not found: " << p->first << endl;
	  abort();
	}
    }

  if (_colsPeriodic.size() % 2 != 0)
    {
	cerr << "BoundaryManager::BoundaryManager()\n";
	cerr << "problem in periodic component data\n";
	cerr << "even number of colors needed\n";
	abort();
    }
}

/*-------------------------------------------------------*/

ostream& BoundaryManager::print(ostream& s) const
{
  s << "BoundaryManager" << endl;
  s << "ColorsEquation     :\t" << _colsEquation << endl;
  s << "ColorsRightHandSide:\t" << _colsRightHandSide << endl;
  s << "ColorsDirichlet    :\t" << _colsDirichlet << endl;
  s << "ComponentsDirichlet:\t" << _compsDirichlet << endl;
  s << "ColorsPeriodic     :\t" << _colsPeriodic << endl;
  s << "ComponentsPeriodic :\t" << _compsPeriodic << endl;
  return s;
}
}
