#ifndef __visudatacompvector_h
#define __visudatacompvector_h

#include  "visudata.h"
#include  "gascoigne.h"

/*----------------------------------------------*/

class VisuDataCompVector : public VisuData
{
 protected:

  const Gascoigne::GlobalVector* _v;

 public:

  VisuDataCompVector();
  VisuDataCompVector(const Gascoigne::GlobalVector& v);

  void SetGlobalVector(const Gascoigne::GlobalVector* v);

  int    visucomp()     const;
  int    visun()        const;
  double visudata(int i,int c) const;
};

/***************************************************************/

#endif
