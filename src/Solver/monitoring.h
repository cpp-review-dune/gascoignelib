#ifndef __monitoring_h
#define __monitoring_h

#include "nvector.h"

/*---------------------------------------------------------*/

namespace Gascoigne
{
class Monitoring
{
 protected:

  std::vector<DoubleVector > Juh, Je;
  DoubleVector          eta, nnodes, ncells;
  DoubleVector          Ju;
  int                      niter;

 public:

  Monitoring() {}
  ~Monitoring() 
    {
      std::cout << "Monitor" << std::endl;
      std::cout << "----------------------------------" << std::endl;
      std::cout << "-- nn   nc   j   eta   je   eff --" << std::endl;
      std::cout << "----------------------------------" << std::endl;
      for (int i=0; i<Juh.size(); i++)
        {
          std::cout << nnodes[i] << " ";
          std::cout << ncells[i] << " ";
          std::cout << Juh[i] << " ";
          std::cout << eta[i] << " ";
          std::cout << Je[i] << " ";
          std::cout << eff(i) << std::endl;
        }
      std::cout << "----------------------------------" << std::endl;
    }
  void BasicInit(const DoubleVector& ju) 
    { 
      Ju = ju; 
    }
  void SetMeshInformation(int iter, int nodes, int cells)
    {
      assert(iter>=1);
      niter = iter;
      Juh.resize(niter);
      Je .resize(niter);
      eta.resize(niter);
      nnodes.resize(niter); nnodes[iter-1] = nodes;
      ncells.resize(niter); ncells[iter-1] = cells;
    }
  void SetSolutionInformation(int iter, DoubleVector juh, double et)
    {
      assert(niter==iter);
      Juh[iter-1] = juh;
      Je [iter-1] = Ju;
      Je [iter-1].add(-1,juh);
      eta[iter-1] = et;
    }
  double eff(int i) const { if(Je[i].size()==0) return -1.; else return eta[i]/Je[i][0];}
  void output() const
    {
      int i = ncells.size()-1;
      std::cout << "## " << nnodes[i] << " ";
      std::cout << ncells[i] << " ";
      std::cout << Juh[i] << " ";
      std::cout << eta[i] << " ";
      std::cout << Je[i] << " ";
      std::cout << eff(i) << std::endl;
    }
};
}

/*---------------------------------------------------------*/

#endif
