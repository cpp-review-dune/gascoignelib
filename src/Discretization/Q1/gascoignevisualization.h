#ifndef  __GascoigneVisualization_h
#define  __GascoigneVisualization_h

#include  "gascoigne.h"
#include  "visualization.h"
#include  "visudatacompvector.h"

/*-----------------------------------------*/

namespace Gascoigne
{
class ComponentInformation;

class GascoigneVisualization : public Visualization
{
protected:

  const GlobalVector* _v;

  VisuDataInfo        VDI;
  VisuDataCompVector  VD;

public:

  GascoigneVisualization() : Visualization(), _v(NULL) {}
  ~GascoigneVisualization() {}

  void AddVector(const ComponentInformation* CI, const GlobalVector* v);
  void AddVector(const GlobalVector* v);
};
}

#endif
