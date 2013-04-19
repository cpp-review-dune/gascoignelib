#include "triacontainer.h"
#include "meshhierarchy.h"
#include "stopwatch.h"
#include "meshagent.h"
#include "hierarchicalmesh2d.h"
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
  MeshAgent MA;
  
  ParamFile pf("old.param");

  HierarchicalMesh2d HM;
  HM.BasicInit(&pf);

  nvector<int> coa;
  
  StopWatch S;
  
  int oldcells=0;
  int newcells=HM.ncells();
  
  cout << "############################## Testing Hierarchical Mesh" 
       << endl << endl;
  
  cout << "step\t" << "cells\t" << "ref\t" << "MB\t" << "B/cell\t"
       << "sec\t" << "sec/Mcell" << endl;
  
  for (int i=0;i<15;++i)
    {
      // refine
      oldcells = newcells;

      S.reset();
      S.start();
      nvector<int> ref;
      for (int jj=0;jj<newcells;++jj)
	if (rand()%200<29)
	  ref.push_back(jj);
      if (i>0)
	HM.patch_refine(ref,coa);
      else
	HM.global_refine(1);
      
      newcells = HM.ncells();
      
      S.stop();


      long int mem = get_mem_usage();
      cout << i << "\t" << newcells << "\t" 
	   << 1.0*newcells/oldcells << "\t"
	   << mem << "\t" << (1024*1024 * mem) / (newcells) << "\t"
	   << S.read() << "\t"
	   << S.read() / (newcells*1.e-6) << endl;
            
      
      
      // S.reset();
      // S.start();
      // MeshHierarchy<__DIM__> MH(TC);
      // MH.ReInit();
      // S.stop();
      // long int hmem = get_mem_usage();
      // cout << "Hierarchy:\t" << MH.nlevels() << " levels" << endl
      // 	   << "\t " << S.read() << "s\t"
      // 	   << 1.e6*S.read()/newcells  << " s/Mcell " << endl
      // 	   << "\t " << hmem-mem << " MB \t"
      // 	   << 1024*1024*(hmem-mem)/newcells << " Byte/cell" << endl;
      // cout << endl;
    }
  
  
  
  return 0;
}
