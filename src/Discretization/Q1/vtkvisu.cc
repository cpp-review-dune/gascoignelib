#include "vtkvisu.h"
#include "visudatanvector.h"
#include "visudatacompvector.h"

namespace Gascoigne
{
/*-------------------------------------------------*/

VtkVisu::~VtkVisu() {}

/*-------------------------------------------------*/

VtkVisu::VtkVisu(const MeshInterface& M, const std::string& name, int iter) : Visualization()
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

void VtkVisu::WriteCellData(const GlobalCellVector& eta)
{
  int ncomp = eta.ncomp();

  VisuDataInfo     VDI(1);
  VDI.Clear();
  //  VisuDataNVector  VD(eta);
  VisuDataCompVector VD(eta);

  VDI.AddScalars(ncomp);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);

  write();
}

/*-------------------------------------------------*/
}
