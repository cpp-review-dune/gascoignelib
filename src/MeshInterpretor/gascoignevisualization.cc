#include  "gascoignevisualization.h"


/*-----------------------------------------*/

void GascoigneVisualization::AddVector(const GlobalVector* v) 
{
  assert(mesh);
  assert(v);
  _v=v;

  int ncomp = v->ncomp();
  assert(ncomp);

  VD .Clear();
  VDI.Clear();

  VD .AddGlobalVector(v);
  VDI.AddScalars(ncomp);

  if (ncomp>2) 
    {
      fixarray<3,int> ff;
      int dim = mesh->dimension();
      if (dim==3)
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = 3;
	}
      else if (dim==2)
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = -1;
	}
      VDI.AddVector("v",ff);
    }

  SetPointData(&VD);
  SetPointDataInfo(&VDI);
}
