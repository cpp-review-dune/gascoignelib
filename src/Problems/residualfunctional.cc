#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"
#include  "constantrighthandside.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
ResidualFunctional::ResidualFunctional() : _DD(NULL) 
{ 
  _comp = 0;
  _scale = 1.;
}

/*-----------------------------------------*/

ResidualFunctional::~ResidualFunctional()
{
  if(_DD!=NULL) {delete _DD; _DD=NULL;}
}

/*-----------------------------------------*/
}
