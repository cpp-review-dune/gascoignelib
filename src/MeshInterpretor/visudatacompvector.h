#ifndef __visudatacompvector_h
#define __visudatacompvector_h

#include  "visudata.h"
#include  "gascoigne.h"

using namespace std;
using namespace Gascoigne;

/*----------------------------------------------*/

// alle ncomps gleich grosss !!!!!!!!!

class VisuDataCompVector : public VisuData
{
 protected:

  vector<const GlobalVector*> vvp;
  vector<int>           isizes;

  pair<int,int> GetIndex(int c) const;


 public:

  virtual ~VisuDataCompVector(){}
  VisuDataCompVector();
  VisuDataCompVector(const GlobalVector& v);

  void Clear() {
    vvp.clear();
    isizes.clear();
  }
  void AddGlobalVector(const GlobalVector* v);
  virtual int    visucomp()     const;
  int    visun()        const;
  virtual double visudata(int i,int c) const;
};

/***************************************************************/

#endif
