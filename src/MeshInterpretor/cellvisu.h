#ifndef __cellvisualization_h
#define __cellvisualization_h

#include "visualization.h"

class CellVisualization : public Visualization
{
 public:

  CellVisualization(const MeshInterface& LM, const nvector<double>& u, 
		    const std::string& name, int i);
  CellVisualization(const MeshInterface& LM, const nvector<int>& u, 
		    const std::string& name, int i);
};

#endif
