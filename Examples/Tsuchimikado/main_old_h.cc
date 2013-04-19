#include "triacontainer.h"
#include "meshhierarchy.h"
#include "stopwatch.h"
#include "meshagent.h"
#include "hierarchicalmesh2d.h"
#include <iostream>
#include <sstream>
#include "gascoignemultigridmesh.h"
#include "gascoignemeshconstructor.h"
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
  
  cout << "############################## Testing old Multigrid Mesh" 
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
         

      S.reset();
      S.start();
      GascoigneMultiGridMesh GMG;
      GMG.ReInit(2,HM.nlevels()-HM.patchdepth());
      GascoigneMeshConstructor MGM(&HM,&GMG);
      MGM.BasicInit();
      S.stop();
      long int hmem = get_mem_usage() - mem;

      cout << i << "\t" << newcells << "\t" 
	   << 1.0*newcells/oldcells << "\t"
	   << hmem << "\t" << (1024*1024 * hmem) / (newcells) << "\t"
	   << S.read() << "\t"
	   << S.read() / (newcells*1.e-6) << endl;

    }
  
  
  
  return 0;
}
