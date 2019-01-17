#include "cgdisc.xx"

#include "finiteelement.h"
#include "elementintegrator.h"
#include "elementlpsintegrator.h"
#include "integrationformula.h"
#include "transformation2d.h"
#include "transformation3d.h"
#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq22d.h"
#include "baseq23d.h"

#include "baseq1patch.h"
#include "baseq13dpatch.h"

#include "patchintegrationformula.h"
#include "finiteelement.xx"
#include "gascoignevisualization.h"
namespace Gascoigne
{
  // Visualization
  template <int DIM, int DEGREE, class FINITEELEMENT, class INTEGRATOR>
  void CGDisc<DIM,DEGREE,FINITEELEMENT,INTEGRATOR>::VisuVtk(const ComponentInformation* CI, const ParamFile& pf,
							    const std::string &name, const GlobalVector& u, int i) const
  {
    HNAverage(const_cast<GlobalVector &>(u));

    GascoigneVisualization Visu;

    Visu.SetMesh(GetDofHandler());
    
    if (CI)
      {
	Visu.AddPointVector(CI, &u);
      }
    else
      {
	Visu.AddPointVector(&u);
      }
    
    Visu.read_parameters(&pf);
    Visu.set_name(name);
    Visu.step(i);
    Visu.write();
    
    HNZero(const_cast<GlobalVector &>(u));
  }
    




  
  template class CGDiscQ12d;
  template class CGDiscQ22d;
  template class CGDiscQ13d;
  template class CGDiscQ23d;

  template class Transformation3d<BaseQ13dPatch>;
  template class FiniteElement<3, 2, Transformation3d<BaseQ13dPatch>, BaseQ13dPatch>;

  template class CGDiscQ13dLps;
  template class CGDiscQ23dLps;

  template class CGDiscQ12dLps;
  template class CGDiscQ22dLps;

}

