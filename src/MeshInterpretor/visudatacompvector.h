#ifndef __visudatacompvector_h
#define __visudatacompvector_h

#include  "visudata.h"
#include  "gascoigne.h"

using namespace std;
using namespace Gascoigne;

/*----------------------------------------------*/

class VisuDataCompVector : public VisuData
{
 protected:

  const GlobalVector* _v;

 public:

  VisuDataCompVector();
  VisuDataCompVector(const GlobalVector& v);

  void SetGlobalVector(const GlobalVector* v);

  virtual int    visucomp()     const;
  int    visun()        const;
  virtual double visudata(int i,int c) const;
};

/***************************************************************/

#endif
