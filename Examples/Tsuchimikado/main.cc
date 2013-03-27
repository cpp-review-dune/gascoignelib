#include "triacontainer.h"
#include "meshhierarchy.h"
#include "stopwatch.h"


#include <iostream>
#include <sstream>

using namespace std;
using namespace Tsuchimikado;
using namespace Gascoigne;


#define __DIM__ 2

#define __REF__ (__DIM__*__DIM__-__DIM__+1)


long int  get_mem_usage ()
{
  struct rusage usage;
  
  
  int who = RUSAGE_SELF;
  getrusage(who, &usage);

  return usage.ru_maxrss/1024;
}


int main(int argc, char** argv)
{
  TriaContainer<__DIM__> TC;//(0);
  
  stringstream meshname;
  meshname << "o" << __DIM__ << "d.inp";
  
  string mn = meshname.str();
  if (argc>=2)
    mn = argv[1];

  TC.read_inp(mn);

  StopWatch S;
  
  int oldcells=0;
  int newcells=TC.ncells();
  
  for (int i=0;i<15;++i)
    {
      // refine
      oldcells = newcells;
      S.reset();
      S.start();
      TC.clear_refine_flags();
      for (int c=0;c<TC.ncells();++c)
	if (rand()%2<1)
	  TC.set_refine_flag(c,3);
      TC.refine_cells();
      newcells = TC.ncells();
      S.stop();
      long int mem = get_mem_usage();
      cout << "Step " << i << " mesh with " << newcells << " total cells" << endl;
      cout << "Refine from \t" << oldcells << "\t -> \t" << newcells
	   << "\t (" << newcells-oldcells << ")" << endl;
      cout << "\t " << S.read() << "s\t "
	   << 1.e6*S.read()/(newcells-oldcells)  << " s/Mcell" << endl
	   << "\t" << mem << " MB \t"
	   << 1024*1024*mem/newcells << " Byte/cell" << endl;
      
      
      
      S.reset();
      S.start();
      MeshHierarchy<__DIM__> MH(TC);
      MH.ReInit();
      S.stop();
      long int hmem = get_mem_usage();
      cout << "Hierarchy:\t" << MH.nlevels() << " levels" << endl
	   << "\t " << S.read() << "s\t"
	   << 1.e6*S.read()/newcells  << " s/Mcell " << endl
	   << "\t " << hmem-mem << " MB \t"
	   << 1024*1024*(hmem-mem)/newcells << " Byte/cell" << endl;
      cout << endl;
    }
  
  
  
  return 0;
}
