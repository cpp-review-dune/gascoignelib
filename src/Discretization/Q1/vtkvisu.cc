#include "vtkvisu.h"
#include "visudatanvector.h"

namespace Gascoigne
{
/*-------------------------------------------------*/

VtkVisu::~VtkVisu() {}

/*-------------------------------------------------*/

VtkVisu::VtkVisu(const MeshInterface& M, std::string name, int iter) : Visualization()
{
  format("vtk");
  set_name(name);
  step(iter);

  SetMesh(&M);
}

/*-------------------------------------------------*/

void VtkVisu::WriteNodeData(const DoubleVector& eta)
{
  VisuDataInfo     VDI(1);
  VisuDataNVector  VD(eta);
  
  SetPointData(&VD);
  SetPointDataInfo(&VDI);

  write();
}

/*-------------------------------------------------*/

void VtkVisu::WriteCellData(const DoubleVector& eta)
{
  VisuDataInfo     VDI(1);
  VisuDataNVector  VD(eta);
  
  SetCellData(&VD);
  SetCellDataInfo(&VDI);

  write();
}

/*-------------------------------------------------*/
}
