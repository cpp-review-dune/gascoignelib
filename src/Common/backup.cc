#include "backup.h"
#include <cassert>
#include <fstream>

using namespace std;
using namespace Gascoigne;

/********************************************************************/

ReadBackUp::ReadBackUp(const string& name, int& size, int& comp)
{
  ifstream file;
  file.open(name.c_str());
  
  assert(file);

  file >> size >> comp;
}

/********************************************************************/

ReadBackUp::ReadBackUp(GlobalVector& u, const string& name)
{
  ifstream file;
  file.open(name.c_str());
  
  assert(file);

  int size, comp;

  file >> size;
  file >> comp;

  //  cout << "BackUp   : reading " << name << ", ";
  //  cout << comp<<" components, "<< size <<" nodes" <<endl;

  if (u.n()!=size)
    {
      cout << "Incompatibility u.n() size " << u.n() << " " << size << endl;
    }
  assert(u.n()==size);

  int v = GascoigneMath::max_int(u.ncomp(),comp);
  
  double d;
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<v; c++)  
	{
	  file >> u(i,c);
	}
      for (int c=v; c<comp ;c++)  { file >> d;}
    }
  string test;
  file >> test;

  assert(test=="BackUpEnd");
}

/********************************************************************/

ReadBackUpResize::ReadBackUpResize(GlobalVector& u, const string& name)
{
  ifstream file;
  file.open(name.c_str());
  
  assert(file);

  int size, comp;

  file >> size;
  file >> comp;

  cout << "BackUp   : reading " << name << ", ";
  cout << comp<<" components, "<< size <<" nodes" <<endl;


  u.ReInit(comp,size);

  if (u.n()!=size)
    {
      cout << "Incompatibility u.n() size " << u.n() << " " << size << endl;
    }
  assert(u.n()==size);

  int v = GascoigneMath::max_int(u.ncomp(),comp);
  
  double d;
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<v; c++)  
	{
	  file >> u(i,c);
	}
      for (int c=v; c<comp ;c++)  { file >> d;}
    }
  string test;
  file >> test;

  assert(test=="BackUpEnd");
}

/********************************************************************/

WriteBackUp::WriteBackUp(const GlobalVector& u, const string& bname)
{
  string name = bname + ".bup";

  ofstream file;
  file.open(name.c_str());
  file.setf(ios::scientific,ios::floatfield);
  
  assert(file);
  file << u.n() << " " << u.ncomp() << endl;

  file.precision(10);
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<u.ncomp(); c++)  
	{
	  file << u(i,c) << " ";
	}
      file << endl;
    }
  file << "BackUpEnd" << endl;
  file.close();
}

/********************************************************************/

WriteBackUpBinary::WriteBackUpBinary(const GlobalVector& u, const string& bname)
{
  string name(bname);
  name += ".bup";
  
  ofstream file;
  file.open(name.c_str());
  
  if(!file)
    {
      cerr << "BackUp: writing error" << endl;
      exit(10);
    }

  u.BinWrite(file);

  file << "BackUpEnd" << endl;
  file.close();
}

/********************************************************************/

ReadBackUpBinary::ReadBackUpBinary(GlobalVector& u, const string& bname)
{
  string name(bname);
  name += ".bup";

  ifstream file;
  file.open(name.c_str());
  
  assert(file);

  u.BinRead(file);

  cout << "BackUp   : reading " << name << ", ";
  cout << u.ncomp() <<" components, "<< u.n() <<" nodes " << endl;

  string test;
  file >> test;
  assert(test=="BackUpEnd");
}

/********************************************************************/

