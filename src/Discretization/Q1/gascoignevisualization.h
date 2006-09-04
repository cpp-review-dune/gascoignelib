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

  void AddVector(const GlobalVector* v);
  void AddVector(const ComponentInformation* CI, const GlobalVector* v);

public:

  GascoigneVisualization() : Visualization(), _v(NULL) {}
  ~GascoigneVisualization() {}

  void AddPointVector(const ComponentInformation* CI, const GlobalVector* v);
  void AddPointVector(const GlobalVector* v);
  void AddCellVector(const ComponentInformation* CI, const GlobalVector* v);
  void AddCellVector(const GlobalVector* v);
};
}

#endif
