#ifndef  __GascoigneVisualization_h
#define  __GascoigneVisualization_h


#include  "gascoigne.h"
#include  "visualization.h"
#include  "visudatacompvector.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class GascoigneVisualization : public Visualization
{
private:

  const GlobalVector* _v;

  VisuDataInfo        VDI;
  VisuDataCompVector  VD;

public:


  GascoigneVisualization() : Visualization(), _v(NULL) {}

  void AddVector(const GlobalVector* v);


};
}

#endif
