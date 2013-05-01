#include "triacontainer.h"
#include "meshhierarchy.h"
#include "stopwatch.h"

#include <iostream>
#include <sstream>

#include "continuouscelldofs.h"

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
  meshname << "m" << __DIM__ << "d.inp";
  
  string mn = meshname.str();
  if (argc>=2)
    mn = argv[1];

  TC.read_inp(mn);
  TC.global_refine();
  
  StopWatch S;
  
  int oldcells=0;
  int newcells=TC.ncells();
  cout << "############################## Testing new MultiGridMesh" 
       << endl << endl;
  
  cout << "step\t" << "cells\t" << "ref\t" <<  "MB\t" << "B/cell\t"
       << "sec\t" << "sec/Mcell" << endl;

  for (int i=0;i<15;++i)
    {
      // refine
      oldcells = newcells;
      S.reset();
      S.start();
      TC.clear_refine_flags();
      if (i>0)
	{
	  for (int c=0;c<TC.ncells();++c)
	    if (rand()%2<1)
	      TC.set_refine_flag(c,3);
	  TC.refine_cells();
	}
      else
	TC.global_refine(1);
      
      
      newcells = TC.ncells();
      
      MeshHierarchy<__DIM__>  MH(TC);
      MH.ReInit();
      S.stop();


      cout << MH.GetMeshLevel(0).size() << "\t";
      
      ContinuousCellDofs<__DIM__,1> CCD;
      CCD.BasicInit(&TC);
      CCD.ReInit(MH.GetMeshLevel(0));
      


      long int hmem = get_mem_usage();
      

      // cout << i << "\t" << newcells << "\t" 
      // 	   << 1.0*newcells/oldcells << "\t"
      // 	   << hmem << "\t" << (1024*1024 * hmem) / (newcells) << "\t"
      // 	   << S.read() << "\t"
      // 	   << S.read() / (newcells*1.e-6) << endl;
    }
  
  
  
  return 0;
}
