#include "backup.h"
#include <cassert>
#include <fstream>

using namespace std;

namespace Gascoigne
{
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
  
  if(!file){
    cerr << "backup file '"<< name << "' not found" << endl;
    assert(file);
    //abort();
  }

  int size, comp;

  file >> size;
  file >> comp;

  //  cout << "BackUp   : reading " << name << ", ";
  //  cout << comp<<" components, "<< size <<" nodes" <<endl;

  if (u.n()!=size)
    {
      cout << "Incompatibility u.n() size u.n()=" << u.n() << " file.n()=" << size << endl;
    }
  assert(u.n()==size);

  int v = max_int(u.ncomp(),comp);
  
  double d;
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<v; c++)  
        {
          double val = 0.;
          file >> val;
          u(i,c) += val;
        }
      for (int c=v; c<comp ;c++)  { file >> d;}
    }
  string test;
  file >> test;

  if(test!="BackUpEnd")
    {
      cout << "error in " <<__FILE__ << ":" << __LINE__ << " : error, test=='"<<test<<"' should be =='BackUpEnd'"<<endl;
      if( u.ncomp()!=comp)
        {
          cout << "probably, because: expected comp nr="<< u.ncomp() <<" NOT EQUAL the bup-file supplied comp nr="<<comp<<endl;
        }
      //assert(test=="BackUpEnd");
      abort();
    }
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

  int v = max_int(u.ncomp(),comp);
  
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

  file.precision(16);
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
}
