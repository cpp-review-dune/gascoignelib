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

void FunctionalManager::Print(std::ostream& os) const
{
  os << "FunctionalManager: " << FF.size() << std::endl;
  for(int i=0;i<FF.size();i++)
    {
      os << "functional " << FF[i]->GetName() << std::endl;
    }
}

/*-----------------------------------------*/

void FunctionalManager::ConstructSet(const std::string& paramfile, const Equation& EQ)
{
  // erstmal names einlesen

  DataFormatHandler DF;
  DF.insert("names",&names);

  FileScanner FS0(DF);
  FS0.NoComplain();
  FS0.readfile(paramfile,"Functional");
  
  int n = names.size();

  // syntax fuer functional: 4 std::strings
  //         name type parameterlist exactvalue grid


  std::vector<std::vector<std::string> > functional(n,std::vector<std::string>(5));
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
      Construct(i, functional[i],EQ);
//       FF[i]->GetName() = names[i];
    }
  gnames.clear();
  for (int i=0; i<n; i++)
    {
      const std::string& grid = functional[i][3];
      if (grid=="yes") gnames.push_back(names[i]);
    }
//   std::cerr << "\tAllFunctional\n";
//   for (int i=0; i<n; i++)
//     {
//       std::cerr << "\t\t" << names[i] << "\t" << FF[i]->GetName()<<std::endl;
//     }
//   std::cerr << "\tGridFunctionals\n" << gnames << std::endl << std::endl;
}

/*-----------------------------------------*/

void FunctionalManager::Construct
(int i, const std::vector<std::string>& functional, const Equation& EQ)
{
  const std::string& type     = functional[0];
  const std::string& params   = functional[1];
  const std::string& exactval = functional[2];
  const std::string& grid     = functional[3];

  FF[i] = ConstructFunctional(type,params,EQ);

  double exact = 0.;
  if (exactval!="unknown") 
    {
      exact = atof(exactval.c_str());
      FF[i]->ExactValue() = exact;
    }
}

/*-----------------------------------------*/

const Functional* FunctionalManager::GetFunctional  (const std::string& name) const 
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
(const std::string& name, const std::string& param, const Equation& EQ)
{
  std::vector<std::string> args = StringSplit(param.c_str(),'_');

  if(name=="zero")
    {
      return new ZeroFunctional;
    }
  if(name=="residual")
    {
      return new ResidualFunctional(EQ, args);
    }
  if(name=="point")
    {
      return new PointFunctional(EQ,args);
    }
  if(name=="domainmean")
    {
      return new DomainMeanFunctional(EQ,args);
    }
  if(name=="constantboundary")
    {
      return new ConstantBoundaryFunctional(EQ,args);
    }
  if(name=="nusselt")
    {
      assert(0);
    }
  std::cerr << "FunctionalManager::ConstructFunctional()\n";
  std::cerr << "unknown functional name: " << name << std::endl;
  abort();
}
