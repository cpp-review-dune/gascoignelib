#include  "functionalmanager.h"
#include  "filescanner.h"
#include  "stlio.h"
#include  "boundaryfunctional.h"
#include  "domainmeanfunctional.h"
#include  "residualfunctional.h"
#include  "pointfunctional.h"
#include  "zerofunctional.h"
#include  "constantboundaryfunctional.h"
#include  "stringutil.h"


using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

FunctionalManager::FunctionalManager()
{
}

/*-----------------------------------------*/

FunctionalManager::~FunctionalManager() 
{
  for(int i=0;i<FF.size();i++)
    {
      if(FF[i])  delete FF[i];  FF[i]=NULL;
    }
}

/*-----------------------------------------*/

void FunctionalManager::Print(ostream& os) const
{
  os << "FunctionalManager: " << FF.size() << endl;
  for(int i=0;i<FF.size();i++)
    {
      os << "functional " << FF[i]->GetName() << endl;
    }
}

/*-----------------------------------------*/

void FunctionalManager::ConstructSet(const ParamFile* paramfile)
{
  // erstmal names einlesen

  DataFormatHandler DF;
  DF.insert("names",&names);

  FileScanner FS0(DF);
  FS0.NoComplain();
  FS0.readfile(paramfile,"Functional");
  
  int n = names.size();

  // syntax fuer functional: 4 strings
  //         name type parameterlist exactvalue grid


  vector<vector<string> > functional(n,vector<string>(5));
  DataFormatHandler DFH;

  DFH.insert("names",&names);
  for (int i=0; i<n; i++)
    {
      DFH.insert(names[i],&functional[i]);
    }
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Functional");
  
  FF.reserve(n);
  FF.resize(n);
  for (int i=0; i<n; i++)
    {
      Construct(i, functional[i]);
//       FF[i]->GetName() = names[i];
    }
  gnames.clear();
  for (int i=0; i<n; i++)
    {
      const string& grid = functional[i][3];
      if (grid=="yes") gnames.push_back(names[i]);
    }
//   cerr << "\tAllFunctional\n";
//   for (int i=0; i<n; i++)
//     {
//       cerr << "\t\t" << names[i] << "\t" << FF[i]->GetName()<<endl;
//     }
//   cerr << "\tGridFunctionals\n" << gnames << endl << endl;
}

/*-----------------------------------------*/

void FunctionalManager::Construct
(int i, const vector<string>& functional)
{
  const string& type     = functional[0];
  const string& params   = functional[1];
  const string& exactval = functional[2];
  const string& grid     = functional[3];

  FF[i] = ConstructFunctional(type,params);

  double exact = 0.;
  if (exactval!="unknown") 
    {
      exact = atof(exactval.c_str());
      FF[i]->ExactValue() = exact;
    }
}

/*-----------------------------------------*/

const Functional* FunctionalManager::GetFunctional  (const string& name) const 
{
  for(int i=0;i<names.size();i++)
    {
      if(names[i]==name)    
	{
	  return FF[i];
	}
    }
  return 0;
}

/*-----------------------------------------*/

Functional* FunctionalManager::ConstructFunctional
(const string& name, const string& param)
{
  vector<string> args = StringSplit(param.c_str(),'_');

  if(name=="zero")
    {
      return new ZeroFunctional;
    }
  if(name=="residual")
    {
      return new ResidualFunctional(args);
    }
  if(name=="point")
    {
      return new PointFunctional(args);
    }
  if(name=="domainmean")
    {
      return new DomainMeanFunctional(args);
    }
  if(name=="constantboundary")
    {
      return new ConstantBoundaryFunctional(args);
    }
  if(name=="nusselt")
    {
      assert(0);
    }
  cerr << "FunctionalManager::ConstructFunctional()\n";
  cerr << "unknown functional name: " << name << endl;
  abort();
}
