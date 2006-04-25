#ifndef __vtkvisu_h
#define __vtkvisu_h

#include  "visualization.h"

/*--------------------------------------------------*/

namespace Gascoigne
{

class VtkVisu : public Visualization
{
 public:

  VtkVisu(const MeshInterface& M, const std::string& name, int iter);
  ~VtkVisu();

  void WriteNodeData(const DoubleVector& eta);
  void WriteCellData(const GlobalVector& eta);
};

/*--------------------------------------------------*/
}

#endif
