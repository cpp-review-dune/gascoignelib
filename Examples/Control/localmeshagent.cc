#include  "localmeshagent.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void LocalMeshAgent::BasicInit(const ParamFile* paramfile)
{
  string gridname("Results/forward.00000.gup");
  int dimension=0;
  int prerefine=0;

  MeshAgent::SetDefaultValues(dimension, gridname, prerefine);
  MeshAgent::ReadParamFile(paramfile);

  if (GetDimension()==2)
    {
      HMP = new HierarchicalMesh2d;
    }
  else if (GetDimension()==3)
    {
      HMP = new HierarchicalMesh3d;
    }
  else
    {
      cout << "dimension of Mesh ? " << GetDimension() << endl;
    }

  int patchdepth=1;
  int epatcher=1;
  HMP->SetParameters(gridname,patchdepth,epatcher);
  assert(HMP);

  GMG = NewMultiGridMesh();

  ReInit();
}
