#include "backup.h"
#include <cassert>
#include <fstream>

/********************************************************************/

ReadBackUp::ReadBackUp(const std::string& name, int& size, int& comp)
{
  std::ifstream file;
  file.open(name.c_str());
  
  assert(file);

  file >> size >> comp;
}

/********************************************************************/

ReadBackUp::ReadBackUp(GlobalVector& u, const std::string& name)
{
  std::ifstream file;
  file.open(name.c_str());
  
  assert(file);

  int size, comp;

  file >> size;
  file >> comp;

  std::cout << "BackUp   : reading " << name << ", ";
  std::cout << comp<<" components, "<< size <<" nodes" <<std::endl;

  if (u.n()!=size)
    {
      std::cout << "Incompatibility u.n() size " << u.n() << " " << size << std::endl;
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
  std::string test;
  file >> test;

  assert(test=="BackUpEnd");
}

/********************************************************************/

ReadBackUpResize::ReadBackUpResize(GlobalVector& u, const std::string& name)
{
  std::ifstream file;
  file.open(name.c_str());
  
  assert(file);

  int size, comp;

  file >> size;
  file >> comp;

  std::cout << "BackUp   : reading " << name << ", ";
  std::cout << comp<<" components, "<< size <<" nodes" <<std::endl;


  u.ReInit(comp,size);

  if (u.n()!=size)
    {
      std::cout << "Incompatibility u.n() size " << u.n() << " " << size << std::endl;
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
  std::string test;
  file >> test;

  assert(test=="BackUpEnd");
}

/********************************************************************/

WriteBackUp::WriteBackUp(const GlobalVector& u, const std::string& bname)
{
  std::string name = bname + ".bup";

  std::ofstream file;
  file.open(name.c_str());
  file.setf(std::ios::scientific,std::ios::floatfield);
  
  assert(file);
  file << u.n() << " " << u.ncomp() << std::endl;

  file.precision(10);
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<u.ncomp(); c++)  
	{
	  file << u(i,c) << " ";
	}
      file << std::endl;
    }
  file << "BackUpEnd" << std::endl;
  file.close();
}

/********************************************************************/

WriteBackUpBinary::WriteBackUpBinary(const GlobalVector& u, const std::string& bname)
{
  std::string name(bname);
  name += ".bup";
  
  std::ofstream file;
  file.open(name.c_str());
  
  if(!file)
    {
      std::cerr << "BackUp: writing error" << std::endl;
      exit(10);
    }

  u.BinWrite(file);

  file << "BackUpEnd" << std::endl;
  file.close();
}

/********************************************************************/

ReadBackUpBinary::ReadBackUpBinary(GlobalVector& u, const std::string& bname)
{
  std::string name(bname);
  name += ".bup";

  std::ifstream file;
  file.open(name.c_str());
  
  assert(file);

  int size, comp;

  u.BinRead(file);

  std::cout << "BackUp   : reading " << name << ", ";
  std::cout << u.ncomp() <<" components, "<< u.n() <<" nodes " << std::endl;

  std::string test;
  file >> test;
  assert(test=="BackUpEnd");
}

/********************************************************************/

