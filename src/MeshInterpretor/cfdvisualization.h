#ifndef __cfdvisualization_h
#define __cfdvisualization_h

#include  "visualization.h"
#include  "visudatacompvector.h"

/*-------------------------------------------------------------------------*/

class CfdVisualization : public Visualization
{
  VisuDataInfo        VDI;
  VisuDataCompVector  VD;

public:

  CfdVisualization(const Gascoigne::GlobalVector& u, const MeshInterface& M) : 
    Visualization(), VDI(u.ncomp()), VD(u)
    {
      SetPointData(&VD);
      SetMesh(M);

      if (u.ncomp()>2) 
	{
	  fixarray<3,int> ff;
	  int dim = M.dimension();
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
      SetPointDataInfo(&VDI);
    }
};

/*-------------------------------------------------------------------------*/

#endif
