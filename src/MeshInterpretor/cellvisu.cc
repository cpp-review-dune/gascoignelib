#include "cellvisu.h"
#include "visudatanvector.h"

using namespace std;

/*************************************************************/

CellVisualization::CellVisualization
(const MeshInterface& M, const nvector<double>& u, const string& name, int i)
{
  VisuDataNVector VD(u);

  VisuDataInfo VDI;

  VDI.AddScalar("eta",0);
  
  SetMesh(M);
  format("vtk");
  set_name(name);
  SetCellData(&VD);
  SetCellDataInfo(&VDI);
  
  step(i);
  write();
}

/*************************************************************/


