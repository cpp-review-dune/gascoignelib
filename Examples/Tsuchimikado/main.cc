#include "triacontainer.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace Tsuchimikado;

#define __DIM__ 3

#define __REF__ (__DIM__*__DIM__-__DIM__+1)


int main(int argc, char** argv)
{
  TriaContainer<__DIM__> TC;//(0);
  
  stringstream meshname;
  meshname << "m" << __DIM__ << "d.inp";
  
  string mn = meshname.str();
  if (argc>=2)
    mn = argv[1];

  TC.read_inp(mn);

  for (int i=0;i<10;++i)
    {
      // refine
      TC.clear_refine_flags();
      
      int r=3;
      
      for (int k=0;k<TC.ncells();++k)
	if (rand()%1==0)
	  TC.set_refine_flag(k,7);//rand()%4);
           
      
      TC.refine_cells();


      // print
      cout << i << ":\t" << TC.nvertices() << "\t" << TC.ncells() << endl;
      stringstream on;
      on << "mesh_" << i << ".gp";
      TC.print_gnuplot(on.str(),0); 
    }
  
  
  return 0;
}
