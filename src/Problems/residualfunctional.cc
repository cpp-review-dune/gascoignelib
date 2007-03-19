#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "constantrighthandside.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  ResidualFunctional::ResidualFunctional() : __DD(NULL) 
  { 
    __comps.clear();
    __scales.clear();
    __cols.clear();
  }
  
  /*-----------------------------------------*/
  
  ResidualFunctional::~ResidualFunctional()
  {
    if(__DD!=NULL) {delete __DD; __DD=NULL;}
  }
  
  /*-----------------------------------------*/
}
