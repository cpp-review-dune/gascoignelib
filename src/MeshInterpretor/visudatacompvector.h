#ifndef __visudatacompvector_h
#define __visudatacompvector_h

#include  "visudata.h"
#include  "gascoigne.h"

/*----------------------------------------------*/

namespace Gascoigne
{
class VisuDataCompVector : public VisuData
{
 protected:

  const GlobalVector* _v;

 public:

  VisuDataCompVector();
  VisuDataCompVector(const GlobalVector& v);

  void SetGlobalVector(const GlobalVector* v);

  int    visucomp()     const;
  int    visun()        const;
  double visudata(int i,int c) const;
};
}

/***************************************************************/

#endif
